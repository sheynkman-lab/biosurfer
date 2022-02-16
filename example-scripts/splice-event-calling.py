# %%
from abc import ABC
from enum import Flag, auto
from itertools import chain
from multiprocessing.sharedctypes import Value
from operator import attrgetter
from typing import Iterable, List, Tuple
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from attrs import frozen, field, evolve
from biosurfer.core.constants import Strand
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger, frozendataclass, get_interval_overlap_graph
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.core.models.nonpersistent import Junction, Position
from graph_tool.topology import label_components, is_bipartite, shortest_path
from graph_tool.draw.cairo_draw import graph_draw
from graph_tool import GraphView
from more_itertools import windowed
from tqdm import tqdm


df = pd.read_csv('complex-splice-events.csv', names=['anchor', 'other', 'score'])
# df['gene'] = df['anchor'].str.split('-').str.get(0)

db = Database('gencode')
session = db.get_session()
txs = Transcript.from_names(session, pd.concat([df['anchor'], df['other']]).unique())

# %%
class SpliceEventCategory(Flag):
    INCLUSION = auto()
    EXCLUSION = auto()
    INTRON = auto()
    DONOR = auto()
    ACCEPTOR = auto()
    EXON = auto()
    UNKNOWN = auto()

    IR = INTRON | INCLUSION
    IX = INTRON | EXCLUSION
    DI = DONOR | INCLUSION
    DX = DONOR | EXCLUSION
    AI = ACCEPTOR | INCLUSION
    AX = ACCEPTOR | EXCLUSION
    EI = EXON | INCLUSION
    ES = EXON | EXCLUSION

    def __repr__(self):
        return self.name

# @dataclass(frozen=True)
@frozen(hash=True)
class BasicSpliceEvent:
    category: 'SpliceEventCategory' = field(default=SpliceEventCategory.UNKNOWN)
    anchor_junctions: tuple['Junction'] = field(factory=tuple)
    other_junctions: tuple['Junction'] = field(factory=tuple)

    @classmethod
    def from_junctions(cls, anchor_junctions: Iterable['Junction'], other_junctions: Iterable['Junction']) -> 'BasicSpliceEvent':
        A = len(anchor_junctions)
        O = len(other_junctions)
        N = A + O
        if N == 1:
            if anchor_junctions:
                category = SpliceEventCategory.IR
            elif other_junctions:
                category = SpliceEventCategory.IX
        elif N == 2:
            anchor_junc = anchor_junctions[0]
            other_junc = other_junctions[0]
            if anchor_junc == other_junc:
                raise ValueError(f'Invalid BasicSpliceEvent between {anchor_junc} and {other_junc}')
            elif anchor_junc.acceptor == other_junc.acceptor:
                if anchor_junc.donor < other_junc.donor:
                    category = SpliceEventCategory.DI
                else:
                    category = SpliceEventCategory.DX
            elif anchor_junc.donor == other_junc.donor:
                if anchor_junc.acceptor > other_junc.acceptor:
                    category = SpliceEventCategory.AI
                else:
                    category = SpliceEventCategory.AX
            else:
                raise ValueError(f'Invalid BasicSpliceEvent between {anchor_junc} and {other_junc}')
        elif N == 3:
            if A == 1:
                category = SpliceEventCategory.EI
            elif O == 1:
                category = SpliceEventCategory.ES
            if len(set(chain.from_iterable((j.donor, j.acceptor) for j in anchor_junctions + other_junctions))) != 4:
                raise ValueError(f'Invalid BasicSpliceEvent between {anchor_junctions} and {other_junctions}')
        else:
            raise ValueError(f'BasicSpliceEvent cannot have more than 3 junctions')
        anchor_juncs = tuple(sorted(anchor_junctions, key=attrgetter('donor')))
        other_juncs = tuple(sorted(other_junctions, key=attrgetter('donor')))
        return BasicSpliceEvent(category, anchor_juncs, other_juncs)
    
    def __invert__(self):
        return evolve(
            self,
            anchor_junctions = self.other_junctions,
            other_junctions = self.anchor_junctions,
            category = self.category ^ SpliceEventCategory.INCLUSION ^ SpliceEventCategory.EXCLUSION
        )
    
    @property
    def delta_nt(self):
        if self.category & SpliceEventCategory.DONOR:
            return self.other_junctions[0].donor - self.anchor_junctions[0].donor
        elif self.category & SpliceEventCategory.ACCEPTOR:
            return self.anchor_junctions[0].acceptor - self.other_junctions[0].acceptor
        elif self.category is SpliceEventCategory.IR:
            return self.anchor_junctions[0].length
        elif self.category is SpliceEventCategory.IX:
            return -self.other_junctions[0].length
        elif self.category is SpliceEventCategory.EI:
            raise NotImplementedError
        elif self.category is SpliceEventCategory.ES:
            raise NotImplementedError
        else:
            raise AttributeError(f'Cannot calculate delta_nt for category {self.category}')
        
            

@frozen(hash=True)
class SpliceEvent:
    category: tuple['SpliceEventCategory']
    members: tuple['BasicSpliceEvent'] = field(repr=False)
    anchor_junctions: tuple['Junction']
    other_junctions: tuple['Junction']

    @classmethod
    def from_basic_events(cls, events: Iterable['BasicSpliceEvent']):
        events = tuple(events)
        category = tuple(event.category for event in events)
        anchor_junctions = sorted({junc for event in events for junc in event.anchor_junctions}, key=attrgetter('donor'))
        other_junctions = sorted({junc for event in events for junc in event.other_junctions}, key=attrgetter('donor'))
        return SpliceEvent(category, events, tuple(anchor_junctions), tuple(other_junctions))

# %%
def call_splice_event(comp: 'GraphView', start: int, stop: int, chr: str, strand: 'Strand') -> 'SpliceEvent':
    def interval_to_junc(interval: tuple['Position', 'Position']) -> 'Junction':
        return Junction(*interval)
    
    anchor_intervals = tuple(comp.vp.label[v] for v in comp.vertices() if comp.vp.transcript[v] == 0)
    other_intervals = tuple(comp.vp.label[v] for v in comp.vertices() if comp.vp.transcript[v] == 1)
    anchor_junctions = tuple(map(interval_to_junc, anchor_intervals))
    other_junctions = tuple(map(interval_to_junc, other_intervals))

    basic_events = []
    if any(g.vp.label[v][0] < start or g.vp.label[v][1] > stop for v in comp.vertices()):
        # TODO: detect altTSS, altPA
        return None
    else:
        try:
            # Call basic splice events
            basic_event = BasicSpliceEvent.from_junctions(anchor_junctions, other_junctions)
            basic_events.append(basic_event)
        except ValueError:
            # Call compound splice events
            N = comp.num_vertices()
            if N == 2:
                donor_event_cat = SpliceEventCategory.DI if anchor_junctions[0].donor < other_junctions[0].donor else SpliceEventCategory.DX
                acceptor_event_cat = SpliceEventCategory.AI if anchor_junctions[0].acceptor > other_junctions[0].acceptor else SpliceEventCategory.AX
                basic_events.append(BasicSpliceEvent(donor_event_cat, anchor_junctions, other_junctions))
                basic_events.append(BasicSpliceEvent(acceptor_event_cat, anchor_junctions, other_junctions))
            elif N > 2:
                first_vertex = min(comp.vertices(), key=lambda v: comp.vp.label[v])
                last_vertex = max(comp.vertices(), key=lambda v: comp.vp.label[v])
                path, _ = shortest_path(comp, first_vertex, last_vertex)
                for r in path[1:-1]:
                    r_junc = interval_to_junc(comp.vp.label[r])
                    neighbors = sorted(comp.vp.label[v] for v in r.out_neighbors())
                    for p, q in windowed(neighbors, 2):
                        p_junc = interval_to_junc(p)
                        q_junc = interval_to_junc(q)
                        event = BasicSpliceEvent(SpliceEventCategory.EI, (r_junc,), (p_junc, q_junc))
                        if comp.vp.transcript[r]:
                            event = ~event
                        basic_events.append(event)

    return SpliceEvent.from_basic_events(basic_events) if basic_events else None
    

# %%
for a, o, _ in df.sort_values('anchor').itertuples(index=False):
    anchor: 'Transcript' = txs[a]
    other: 'Transcript' = txs[o]
    anchor_intervals = {(junc.donor, junc.acceptor) for junc in anchor.junctions}
    other_intervals = {(junc.donor, junc.acceptor) for junc in other.junctions}
    common_intervals = anchor_intervals & other_intervals
    diff_intervals = anchor_intervals ^ other_intervals

    g, label_to_vertex = get_interval_overlap_graph(((p, q+1) for p, q in diff_intervals), diff_intervals, label_type='object')
    if not is_bipartite(g):
        warn(f'Overlap graph not bipartite for {a} | {o}')
    g.vp.transcript = g.new_vertex_property('int')
    for v in g.vertices():
        g.vp.transcript[v] = g.vp.label[v] in other_intervals
    # for interval in common_intervals:
    #     v = g.add_vertex()
    #     g.add_edge(v, v)
    #     g.vp.label[v] = interval
    #     g.vp.transcript[v] = 2

    g.vp.comp, hist = label_components(g)
    components = {c: GraphView(g, vfilt=lambda v: g.vp.comp[v] == c) for c in range(len(hist))}
    
    print(f'{anchor.name} | {other.name}')
    # graph_draw(g, vertex_text=g.vp.label.coerce_type('string'), vertex_fill_color=g.vp.transcript, vertex_size=10, edge_color='#dddddd')
    chr = anchor.chromosome.name
    strand = anchor.strand
    start = Position(chr, strand, max(anchor.start, other.start))
    stop = Position(chr, strand, min(anchor.stop, other.stop))
    if start > stop:
        start, stop = stop, start
    events = sorted(
        filter(None, (call_splice_event(comp, start, stop, chr, strand) for comp in components.values())),
        key = lambda e: min(j.donor for j in (e.anchor_junctions + e.other_junctions))
    )
    for event in events:
        if event.category:
            print(event)

# %%
