# %%
import heapq
import warnings
from enum import Enum, auto
from itertools import chain
from operator import attrgetter, methodcaller
from typing import Iterable, List, Tuple
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from attrs import evolve, field, frozen
from biosurfer.core.database import Database
from biosurfer.core.helpers import (ExceptionLogger,
                                    get_interval_overlap_graph)
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.core.models.nonpersistent import Junction, Position
from graph_tool import GraphView
from graph_tool.draw.cairo_draw import graph_draw
from graph_tool.topology import is_bipartite, label_components, shortest_path
from more_itertools import windowed
from tqdm import tqdm

df = pd.read_csv('complex-splice-events.csv', names=['anchor', 'other', 'score'])
# df['gene'] = df['anchor'].str.split('-').str.get(0)

db = Database('gencode')
session = db.get_session()
txs = Transcript.from_names(session, pd.concat([df['anchor'], df['other']]).unique())

# %%
class SpliceEventCategory(Enum):
    IR = auto()
    DI = auto()
    AI = auto()
    EI = auto()
    IX = -IR
    DX = -DI
    AX = -AI
    ES = -EI
    UNKNOWN = 0
    
    def __repr__(self):
        return self.name
    
    def __neg__(self):
        return SpliceEventCategory(-self.value)


INCLUSION = {SpliceEventCategory.IR, SpliceEventCategory.DI, SpliceEventCategory.AI, SpliceEventCategory.EI}
EXCLUSION = {SpliceEventCategory.IX, SpliceEventCategory.DX, SpliceEventCategory.AX, SpliceEventCategory.ES}
INTRON = {SpliceEventCategory.IR, SpliceEventCategory.IX}
DONOR = {SpliceEventCategory.DI, SpliceEventCategory.DX}
ACCEPTOR = {SpliceEventCategory.AI, SpliceEventCategory.AX}
EXON = {SpliceEventCategory.EI, SpliceEventCategory.ES}


@frozen(hash=True)
class BasicSpliceEvent:
    category: 'SpliceEventCategory' = field(default=SpliceEventCategory.UNKNOWN)
    anchor_junctions: tuple['Junction'] = field(factory=tuple)
    other_junctions: tuple['Junction'] = field(factory=tuple)
    
    def __neg__(self):
        return evolve(
            self,
            anchor_junctions = self.other_junctions,
            other_junctions = self.anchor_junctions,
            category = -self.category
        )
    
    @property
    def delta_nt(self):
        if self.category in DONOR:
            return self.other_junctions[0].donor - self.anchor_junctions[0].donor
        elif self.category in ACCEPTOR:
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
def call_splice_event(comp: 'GraphView', start: 'Position', stop: 'Position') -> 'SpliceEvent':
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
        if comp.num_vertices() == 1:
            # call alt. intron event
            intron_event_cat = SpliceEventCategory.IR if anchor_junctions else SpliceEventCategory.IX
            basic_events = [BasicSpliceEvent(intron_event_cat, anchor_junctions, other_junctions)]
        else:
            # check for alt. donor usage
            v0, v1 = heapq.nsmallest(2, comp.vertices(), key=lambda v: comp.vp.label[v])
            junc0 = interval_to_junc(comp.vp.label[v0])
            junc1 = interval_to_junc(comp.vp.label[v1])
            if junc0.donor != junc1.donor:
                anchor_junc, other_junc = (junc0, junc1) if comp.vp.transcript[v0] == 0 else (junc1, junc0)
                donor_event_cat = SpliceEventCategory.DI if anchor_junc.donor < other_junc.donor else SpliceEventCategory.DX
                donor_event = [BasicSpliceEvent(donor_event_cat, (anchor_junc,), (other_junc,))]
            else:
                donor_event = []
            # check for alt. acceptor usage
            vM, vN = heapq.nlargest(2, comp.vertices(), key=lambda v: comp.vp.label[v][::-1])
            juncM = interval_to_junc(comp.vp.label[vM])
            juncN = interval_to_junc(comp.vp.label[vN])
            if juncM.acceptor != juncN.acceptor:
                anchor_junc, other_junc = (juncM, juncN) if comp.vp.transcript[vM] == 0 else (juncN, juncM)
                donor_event_cat = SpliceEventCategory.AI if anchor_junc.acceptor > other_junc.acceptor else SpliceEventCategory.AX
                acceptor_event = [BasicSpliceEvent(donor_event_cat, (anchor_junc,), (other_junc,))]
            else:
                acceptor_event = []
            # list all alt. exon events
            exon_events = []
            path, _ = shortest_path(
                comp,
                source = min(v0, v1, key=methodcaller('out_degree')),
                target = min(vM, vN, key=methodcaller('out_degree'))
            )
            for r in path[1:-1]:
                r_junc = interval_to_junc(comp.vp.label[r])
                neighbors = sorted(comp.vp.label[v] for v in r.out_neighbors())
                for p, q in windowed(neighbors, 2):
                    p_junc = interval_to_junc(p)
                    q_junc = interval_to_junc(q)
                    exon_event = BasicSpliceEvent(SpliceEventCategory.EI, (r_junc,), (p_junc, q_junc))
                    if comp.vp.transcript[r] == 1:
                        exon_event = -exon_event
                    exon_events.append(exon_event)
            basic_events = donor_event + exon_events + acceptor_event

    return SpliceEvent.from_basic_events(basic_events) if basic_events else None
    

# %%
all_events = dict()
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
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        events = sorted(
            filter(None, (call_splice_event(comp, start, stop) for comp in components.values())),
            key = lambda e: min(j.donor for j in (e.anchor_junctions + e.other_junctions))
        )
    all_events[a, o] = events
    for event in events:
        if event.category:
            print(f'\t{event.category}')
            for j in event.anchor_junctions:
                print('\t\t' + str(j))
            print('\t\t---')
            for j in event.other_junctions:
                print('\t\t' + str(j))

# %%
