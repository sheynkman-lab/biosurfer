# %%
from abc import ABC, abstractmethod
import heapq
import warnings
from operator import attrgetter, methodcaller
from typing import Iterable
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
@frozen(hash=True)
class BasicSpliceEvent(ABC):
    is_deletion: bool = False

    def __neg__(self):
        return evolve(self, is_deletion=not self.is_deletion)

    @property
    def is_insertion(self):
        return not self.is_deletion

    @property
    @abstractmethod
    def anchor_junctions(self) -> tuple['Junction']:
        pass

    @property
    @abstractmethod
    def other_junctions(self) -> tuple['Junction']:
        pass

    @property
    @abstractmethod
    def delta_nt(self) -> int:
        pass


@frozen(hash=True)
class IntronSpliceEvent(BasicSpliceEvent):
    junction: 'Junction' = field(default=None)

    @property
    def anchor_junctions(self) -> tuple['Junction']:
        return () if self.is_deletion else (self.junction,)
    
    @property
    def other_junctions(self) -> tuple['Junction']:
        return (self.junction,) if self.is_deletion else ()

    @property
    def delta_nt(self) -> int:
        return (-1 if self.is_deletion else 1) * self.junction.length


@frozen(hash=True)
class DonorSpliceEvent(BasicSpliceEvent):
    upstream_junction: 'Junction' = field(default=None)
    downstream_junction: 'Junction' = field(default=None)
    @downstream_junction.validator
    def _check_junctions(self, attribute, value: 'Junction'):
        if self.upstream_junction.donor >= value.donor:
            raise ValueError(f'{value} has donor upstream of {self.upstream_junction}')
    
    @property
    def anchor_junctions(self) -> tuple['Junction']:
        return (self.downstream_junction,) if self.is_deletion else (self.upstream_junction,)

    @property
    def other_junctions(self) -> tuple['Junction']:
        return (self.upstream_junction,) if self.is_deletion else (self.downstream_junction,)

    @property
    def delta_nt(self):
        return (-1 if self.is_deletion else 1) * (self.downstream_junction.donor - self.upstream_junction.donor)


@frozen(hash=True)
class AcceptorSpliceEvent(BasicSpliceEvent):
    upstream_junction: 'Junction' = field(default=None)
    downstream_junction: 'Junction' = field(default=None)
    @downstream_junction.validator
    def _check_junctions(self, attribute, value: 'Junction'):
        if self.upstream_junction.acceptor >= value.acceptor:
            raise ValueError(f'{value} has acceptor upstream of {self.upstream_junction}')
    
    @property
    def anchor_junctions(self) -> tuple['Junction']:
        return (self.upstream_junction,) if self.is_deletion else (self.downstream_junction,)

    @property
    def other_junctions(self) -> tuple['Junction']:
        return (self.downstream_junction,) if self.is_deletion else (self.upstream_junction,)

    @property
    def delta_nt(self):
        return (-1 if self.is_deletion else 1) * (self.downstream_junction.acceptor - self.upstream_junction.acceptor)


@frozen(hash=True)
class ExonSpliceEvent(BasicSpliceEvent):
    skip_junction: 'Junction' = field(default=None)
    upstream_junction: 'Junction' = field(default=None)
    downstream_junction: 'Junction' = field(default=None)
    @downstream_junction.validator
    def _check_short_junctions(self, attribute, value: 'Junction'):
        if self.upstream_junction & value:
            raise ValueError(f'{self.upstream_junction} and {value} overlap')
    
    @property
    def anchor_junctions(self) -> tuple['Junction']:
        return (self.upstream_junction, self.downstream_junction) if self.is_deletion else (self.skip_junction,)
    
    @property
    def other_junctions(self) -> tuple['Junction']:
        return (self.skip_junction,) if self.is_deletion else (self.upstream_junction, self.downstream_junction)

    @property
    def delta_nt(self) -> int:
        return (-1 if self.is_deletion else 1) * (self.downstream_junction.donor - self.upstream_junction.acceptor - 1)


SPLICE_EVENT_CODES = {
    IntronSpliceEvent: ('I', 'i'),
    DonorSpliceEvent: ('D', 'd'),
    AcceptorSpliceEvent: ('A', 'a'),
    ExonSpliceEvent: ('E', 'e')
}


@frozen(hash=True)
class SpliceEvent:
    code: str
    members: tuple['BasicSpliceEvent'] = field(repr=False)
    anchor_junctions: tuple['Junction']
    other_junctions: tuple['Junction']

    @members.validator
    def _check_members(self, attribute, value):
        if len(value) > 1:
            if any(isinstance(event, IntronSpliceEvent) for event in value):
                raise ValueError('Cannot combine IntronSpliceEvent with other BasicSpliceEvents')
            if any(isinstance(event, DonorSpliceEvent) for event in value[1:]):
                raise ValueError('DonorSpliceEvent must be first')
            if any(isinstance(event, AcceptorSpliceEvent) for event in value[:-1]):
                raise ValueError('AcceptorSpliceEvent must be last')

    @classmethod
    def from_basic_events(cls, events: Iterable['BasicSpliceEvent']):
        events = tuple(events)
        code = ''.join(SPLICE_EVENT_CODES[type(event)][event.is_deletion] for event in events)
        anchor_junctions = sorted({junc for event in events for junc in event.anchor_junctions}, key=attrgetter('donor'))
        other_junctions = sorted({junc for event in events for junc in event.other_junctions}, key=attrgetter('donor'))
        return SpliceEvent(code, events, tuple(anchor_junctions), tuple(other_junctions))

# %%
def call_transcript_event(comp: 'GraphView', start: 'Position', stop: 'Position') -> 'SpliceEvent':
    def by_donor(v):
        junc = comp.vp.label[v]
        return junc.donor, junc.acceptor

    def by_acceptor(v):
        junc = comp.vp.label[v]
        return junc.acceptor, junc.donor

    basic_events = []
    N = comp.num_vertices()
    if any(g.vp.label[v].donor < start or g.vp.label[v].acceptor > stop for v in comp.vertices()):
        # TODO: detect altTSS, altPA
        return None
    else:
        if N == 1:
            # call alt. intron event
            v = next(comp.vertices())
            basic_events = [IntronSpliceEvent(is_deletion=not comp.vp.from_anchor[v], junction=comp.vp.label[v])]
        else:
            v0, v1 = heapq.nsmallest(2, comp.vertices(), key=by_donor)
            vN, vM = heapq.nlargest(2, comp.vertices(), key=by_acceptor)
            # check for alt. donor usage
            junc0 = comp.vp.label[v0]
            junc1 = comp.vp.label[v1]
            if junc0.donor != junc1.donor:
                donor_event = DonorSpliceEvent(
                    is_deletion = comp.vp.from_anchor[v1],
                    upstream_junction = junc0,
                    downstream_junction = junc1
                )
                donor_event = [donor_event]
            else:
                donor_event = []
            # check for alt. acceptor usage
            juncM = comp.vp.label[vM]
            juncN = comp.vp.label[vN]
            if juncM.acceptor != juncN.acceptor:
                acceptor_event = AcceptorSpliceEvent(
                    is_deletion = comp.vp.from_anchor[vM],
                    upstream_junction = juncM,
                    downstream_junction = juncN
                )
                acceptor_event = [acceptor_event]
            else:
                acceptor_event = []
            # call any alt. exon events
            exon_events = []
            path, _ = shortest_path(
                comp,
                source = min(v0, v1, key=methodcaller('out_degree')),
                target = min(vM, vN, key=methodcaller('out_degree'))
            )
            for skip in path[1:-1]:
                neighbors = sorted(skip.out_neighbors(), key=by_donor)
                other_has_skip = not comp.vp.from_anchor[skip]
                for upstream, downstream in windowed(neighbors, 2):
                    exon_event = ExonSpliceEvent(
                        is_deletion = other_has_skip,
                        upstream_junction = comp.vp.label[upstream],
                        downstream_junction = comp.vp.label[downstream],
                        skip_junction = comp.vp.label[skip]
                    )
                    exon_events.append(exon_event)
            basic_events = donor_event + exon_events + acceptor_event
        return SpliceEvent.from_basic_events(basic_events) if basic_events else None
    

# %%
all_events = dict()
for a, o, _ in df.sort_values('anchor').itertuples(index=False):
    anchor: 'Transcript' = txs[a]
    other: 'Transcript' = txs[o]
    anchor_junctions = set(anchor.junctions)
    common_junctions = set(anchor.junctions) & set(other.junctions)
    diff_junctions = set(anchor.junctions) ^ set(other.junctions)

    g, label_to_vertex = get_interval_overlap_graph(((j.donor, j.acceptor+1) for j in diff_junctions), diff_junctions, label_type='object')
    if not is_bipartite(g):
        warn(f'Overlap graph not bipartite for {a} | {o}')
    g.vp.from_anchor = g.new_vertex_property('bool')
    for v in g.vertices():
        g.vp.from_anchor[v] = g.vp.label[v] in anchor_junctions
    # for interval in common_intervals:
    #     v = g.add_vertex()
    #     g.add_edge(v, v)
    #     g.vp.label[v] = interval
    #     g.vp.from_anchor[v] = 2

    g.vp.comp, hist = label_components(g)
    components = {c: GraphView(g, vfilt=lambda v: g.vp.comp[v] == c) for c in range(len(hist))}
    
    print(f'{anchor.name} | {other.name}')
    g.vp.text = g.new_vertex_property('string')
    for v in g.vertices():
        junc = g.vp.label[v]
        g.vp.text[v] = f'({junc.donor.coordinate % 10**5}, {junc.acceptor.coordinate % 10**5})'
    # graph_draw(g, vertex_text=g.vp.text, vertex_fill_color=g.vp.from_anchor.coerce_type('int'), vertex_size=8, edge_color='#dddddd')
    chr = anchor.chromosome.name
    strand = anchor.strand
    start = Position(chr, strand, max(anchor.start, other.start))
    stop = Position(chr, strand, min(anchor.stop, other.stop))
    if start > stop:
        start, stop = stop, start
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        events = sorted(
            filter(None, (call_transcript_event(comp, start, stop) for comp in components.values())),
            key = lambda e: min(j.donor for j in (e.anchor_junctions + e.other_junctions))
        )
    all_events[a, o] = events
    for event in events:
        if event.code:
            print(f'\t{event.code}')
            for j in event.anchor_junctions:
                print('\t\t' + str(j))
            print('\t\t---')
            for j in event.other_junctions:
                print('\t\t' + str(j))

# %%
