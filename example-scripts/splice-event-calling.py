# %%
from collections import deque
from functools import cached_property
import heapq
from itertools import groupby
import os
from time import time
import warnings
from abc import ABC, abstractmethod
from operator import attrgetter, methodcaller
from typing import Any, Iterable, Optional, Sequence, Union
from warnings import warn

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from attrs import define, evolve, field, frozen, validators
from biosurfer.core.constants import TranscriptLevelAlignmentCategory
from biosurfer.core.database import Database
from biosurfer.core.helpers import (BisectDict, ExceptionLogger, Interval,
                                    IntervalTree, get_interval_overlap_graph)
from biosurfer.core.models.biomolecules import Exon, Transcript
from biosurfer.core.models.nonpersistent import Junction, Position
from biosurfer.plots.plotting import IsoformPlot
from graph_tool import GraphView
from graph_tool.topology import is_bipartite, label_components, shortest_path
from more_itertools import first, partition, windowed
from tqdm import tqdm

df = pd.read_csv('complex-splice-events.csv', names=['anchor', 'other', 'score'])
# df['gene'] = df['anchor'].str.split('-').str.get(0)

db = Database('gencode')
session = db.get_session()
txs = Transcript.from_names(session, pd.concat([df['anchor'], df['other']]).unique())

# %%
@frozen(eq=True)
class TranscriptEvent(ABC):
    @abstractmethod
    def __neg__(self) -> 'TranscriptEvent':
        raise NotImplementedError

    @property
    @abstractmethod
    def delta_nt(self) -> int:
        raise NotImplementedError
    
    @property
    @abstractmethod
    def start(self) -> 'Position':
        raise NotImplementedError
    
    @property
    @abstractmethod
    def stop(self) -> 'Position':
        raise NotImplementedError


def sort_events(x: Iterable['TranscriptEvent']):
    return tuple(sorted(x, key=attrgetter('start', 'stop')))


@frozen(eq=True)
class BasicTranscriptEvent(TranscriptEvent):
    is_deletion: bool

    def __neg__(self):
        return evolve(self, is_deletion=not self.is_deletion)

    @property
    def is_insertion(self):
        return not self.is_deletion

    @property
    def delta_nt(self) -> int:
        return (-1 if self.is_deletion else 1) * self.length

    @property
    def length(self) -> int:
        return (self.stop - self.start) + 1


@frozen(eq=True)
class CompoundTranscriptEvent(TranscriptEvent):
    members: tuple['BasicTranscriptEvent', ...] = field(converter=sort_events)

    def __neg__(self):
        return self.from_basic_events(-event for event in self.members)

    @property
    def delta_nt(self) -> int:
        return sum(event.delta_nt for event in self.members)

    @property
    def start(self):
        return self.members[0].start

    @property
    def stop(self):
        return self.members[-1].stop

    @classmethod
    def from_basic_events(cls, basic_events: Iterable['BasicTranscriptEvent']):
        return cls(members=basic_events)


@frozen(eq=True)
class IntronSpliceEvent(BasicTranscriptEvent):
    junction: 'Junction'

    @property
    def anchor_junctions(self):
        return () if self.is_deletion else (self.junction,)
    
    @property
    def other_junctions(self):
        return (self.junction,) if self.is_deletion else ()

    @property
    def start(self) -> 'Position':
        return self.junction.donor
    
    @property
    def stop(self) -> 'Position':
        return self.junction.acceptor


@frozen(eq=True)
class DonorSpliceEvent(BasicTranscriptEvent):
    upstream_junction: 'Junction'
    downstream_junction: 'Junction' = field()
    @downstream_junction.validator
    def _check_junctions(self, attribute, value: 'Junction'):
        if self.upstream_junction.donor >= value.donor:
            raise ValueError(f'{value} has donor upstream of {self.upstream_junction}')
    
    @property
    def anchor_junctions(self):
        return (self.downstream_junction,) if self.is_deletion else (self.upstream_junction,)

    @property
    def other_junctions(self):
        return (self.upstream_junction,) if self.is_deletion else (self.downstream_junction,)

    @property
    def start(self) -> 'Position':
        return self.upstream_junction.donor
    
    @property
    def stop(self) -> 'Position':
        return self.downstream_junction.donor - 1


@frozen(eq=True)
class AcceptorSpliceEvent(BasicTranscriptEvent):
    upstream_junction: 'Junction'
    downstream_junction: 'Junction' = field()
    @downstream_junction.validator
    def _check_junctions(self, attribute, value: 'Junction'):
        if self.upstream_junction.acceptor >= value.acceptor:
            raise ValueError(f'{value} has acceptor upstream of {self.upstream_junction}')
    
    @property
    def anchor_junctions(self):
        return (self.upstream_junction,) if self.is_deletion else (self.downstream_junction,)

    @property
    def other_junctions(self):
        return (self.downstream_junction,) if self.is_deletion else (self.upstream_junction,)

    @property
    def start(self) -> 'Position':
        return self.upstream_junction.acceptor + 1
    
    @property
    def stop(self) -> 'Position':
        return self.downstream_junction.acceptor


@frozen(eq=True)
class ExonSpliceEvent(BasicTranscriptEvent):
    skip_junction: 'Junction'
    upstream_junction: 'Junction'
    downstream_junction: 'Junction' = field()
    @downstream_junction.validator
    def _check_short_junctions(self, attribute, value: 'Junction'):
        if self.upstream_junction & value:
            raise ValueError(f'{self.upstream_junction} and {value} overlap')
    
    @property
    def anchor_junctions(self):
        return (self.upstream_junction, self.downstream_junction) if self.is_deletion else (self.skip_junction,)
    
    @property
    def other_junctions(self):
        return (self.skip_junction,) if self.is_deletion else (self.upstream_junction, self.downstream_junction)

    @property
    def start(self) -> 'Position':
        return self.upstream_junction.acceptor + 1
    
    @property
    def stop(self) -> 'Position':
        return self.downstream_junction.donor - 1


BasicSpliceEvent = Union[IntronSpliceEvent, DonorSpliceEvent, AcceptorSpliceEvent, ExonSpliceEvent]


SPLICE_EVENT_CODES = {
    IntronSpliceEvent: ('I', 'i'),
    DonorSpliceEvent: ('D', 'd'),
    AcceptorSpliceEvent: ('A', 'a'),
    ExonSpliceEvent: ('E', 'e')
}


@frozen(eq=True)
class SpliceEvent(CompoundTranscriptEvent):
    members: tuple['BasicSpliceEvent', ...] = field(converter=sort_events, repr=False)
    code: str = field(default='')
    anchor_junctions: tuple['Junction', ...] = field(factory=tuple)
    other_junctions: tuple['Junction', ...] = field(factory=tuple)

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
        return cls(members=events, code=code, anchor_junctions=tuple(anchor_junctions), other_junctions=tuple(other_junctions))


@frozen(eq=True)
class ExonBypassEvent(BasicTranscriptEvent):
    exon: 'Junction'  # TODO: replace with updated Exon object
    is_partial: bool = False

    @property
    def start(self) -> 'Position':
        return self.exon.donor
    
    @property
    def stop(self) -> 'Position':
        return self.exon.acceptor


@frozen(eq=True)
class TSSEvent(CompoundTranscriptEvent):
    members: tuple['ExonBypassEvent', ...] = field(converter=sort_events)
    
    @members.validator
    def _check_members(self, attribute, value: tuple['ExonBypassEvent', ...]):
        if any(event.is_partial for event in value[:-1]):
            raise ValueError
        if len({event.is_deletion for event in value[:-1]}) > 1:
            raise ValueError


@frozen(eq=True)
class APAEvent(CompoundTranscriptEvent):
    members: tuple['ExonBypassEvent', ...] = field(converter=sort_events)

    @members.validator
    def _check_members(self, attribute, value: tuple['ExonBypassEvent', ...]):
        if any(event.is_partial for event in value[1:]):
            raise ValueError
        if len({event.is_deletion for event in value[1:]}) > 1:
            raise ValueError


def call_splice_event(comp: 'GraphView') -> 'SpliceEvent':
    def by_donor(v):
        junc = comp.vp.label[v]
        return junc.donor, junc.acceptor

    def by_acceptor(v):
        junc = comp.vp.label[v]
        return junc.acceptor, junc.donor

    basic_events = []
    N = comp.num_vertices()
    if N == 1:
        # call alt. intron event
        v = first(comp.vertices())
        if not (comp.vp.overlaps_tss[v] or comp.vp.overlaps_pas[v]):
            basic_events = [IntronSpliceEvent(is_deletion=not comp.vp.from_anchor[v], junction=comp.vp.label[v])]
    else:
        v0, v1 = heapq.nsmallest(2, comp.vertices(), key=by_donor)
        vN, vM = heapq.nlargest(2, comp.vertices(), key=by_acceptor)
        # check for alt. donor usage
        junc0 = comp.vp.label[v0]
        junc1 = comp.vp.label[v1]
        if junc0.donor != junc1.donor and not comp.vp.overlaps_tss[v0]:
            donor_event = DonorSpliceEvent(
                is_deletion = bool(comp.vp.from_anchor[v1]),
                upstream_junction = junc0,
                downstream_junction = junc1
            )
            donor_event = [donor_event]
        else:
            donor_event = []
        # check for alt. acceptor usage
        juncM = comp.vp.label[vM]
        juncN = comp.vp.label[vN]
        if juncM.acceptor != juncN.acceptor and not comp.vp.overlaps_pas[vN]:
            acceptor_event = AcceptorSpliceEvent(
                is_deletion = bool(comp.vp.from_anchor[vM]),
                upstream_junction = juncM,
                downstream_junction = juncN
            )
            acceptor_event = [acceptor_event]
        else:
            acceptor_event = []
        # call any alt. exon events
        exon_events = []
        if N > 2:
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


def call_transcript_events(anchor: 'Transcript', other: 'Transcript'):
    chr = anchor.chromosome.name
    strand = anchor.strand
    anchor_start = Position(chr, strand, anchor.start)
    anchor_stop = Position(chr, strand, anchor.stop)
    if anchor_start > anchor_stop:
        anchor_start, anchor_stop = anchor_stop, anchor_start
    other_start = Position(chr, strand, other.start)
    other_stop = Position(chr, strand, other.stop)
    if other_start > other_stop:
        other_start, other_stop = other_stop, other_start
    downstream_start = max(anchor_start, other_start)
    upstream_stop = min(anchor_stop, other_stop)

    anchor_junctions = set(anchor.junctions)
    diff_junctions = anchor_junctions ^ set(other.junctions)
    diff_junctions = {junc for junc in diff_junctions if downstream_start <= junc.acceptor and junc.donor <= upstream_stop}
    tss_overlap_junction = first((junc for junc in diff_junctions if junc.donor <= downstream_start <= junc.acceptor + 1), None)
    pas_overlap_junction = first((junc for junc in diff_junctions if junc.donor - 1 <= upstream_stop <= junc.acceptor), None)

    g, _ = get_interval_overlap_graph(((j.donor, j.acceptor+1) for j in diff_junctions), diff_junctions, label_type='object')
    if not is_bipartite(g):
        warn(f'Overlap graph not bipartite for {a} | {o}')
    g.vp.from_anchor = g.new_vertex_property('bool')
    g.vp.overlaps_tss = g.new_vertex_property('bool')
    g.vp.overlaps_pas = g.new_vertex_property('bool')
    for v in g.vertices():
        junc = g.vp.label[v]
        g.vp.from_anchor[v] = junc in anchor_junctions
        g.vp.overlaps_tss[v] = junc == tss_overlap_junction
        g.vp.overlaps_pas[v] = junc == pas_overlap_junction

    g.vp.comp, hist = label_components(g)
    components = {c: GraphView(g, vfilt=lambda v: g.vp.comp[v] == c) for c in range(len(hist))}
    
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        splice_events = sorted(
            filter(None, (call_splice_event(comp) for comp in components.values())),
            key = lambda e: min(j.donor for j in (e.anchor_junctions + e.other_junctions))
        )
    
    # TODO: simplify this when Exons are refactored
    def get_exon(exon_obj: 'Exon'):
        return Junction(*sorted((Position(chr, strand, exon_obj.start), Position(chr, strand, exon_obj.stop))))
    
    anchor_exons = [get_exon(exon) for exon in anchor.exons]
    other_exons = [get_exon(exon) for exon in other.exons]
    upstream_exons = sorted((exon for exon in set(anchor_exons) | set(other_exons) if exon.donor < downstream_start), key=attrgetter('donor'))
    downstream_exons = sorted((exon for exon in set(anchor_exons) | set(other_exons) if upstream_stop < exon.acceptor), key=attrgetter('donor'))
    # call TSS event
    if upstream_exons:
        is_deletion = downstream_start == other_start
        bypass_events = [ExonBypassEvent(is_deletion, exon) for exon in upstream_exons]
        if tss_overlap_junction:
            alt_downstream_exon = first(other_exons if is_deletion else anchor_exons)
            last_bypass_event = ExonBypassEvent(
                not is_deletion,
                Junction(
                    downstream_start,
                    min(alt_downstream_exon.acceptor, tss_overlap_junction.acceptor)
                ),
                is_partial = tss_overlap_junction.acceptor < alt_downstream_exon.acceptor
            )
            bypass_events.append(last_bypass_event)
        elif any(exon.donor < downstream_start <= exon.acceptor for exon in (anchor_exons if is_deletion else other_exons)):
            upstream_exons[-1] = evolve(upstream_exons[-1], acceptor=downstream_start-1)
            bypass_events[-1] = evolve(bypass_events[-1], exon=upstream_exons[-1], is_partial=True)
        tss_event = TSSEvent.from_basic_events(bypass_events)
    else:
        tss_event = None
    # call APA event
    if downstream_exons:
        is_deletion = upstream_stop == other_stop
        bypass_events = [ExonBypassEvent(is_deletion, exon) for exon in downstream_exons]
        if pas_overlap_junction:
            alt_upstream_exon = first(other_exons if is_deletion else anchor_exons)
            first_bypass_event = ExonBypassEvent(
                not is_deletion,
                Junction(
                    max(alt_upstream_exon.donor, pas_overlap_junction.donor),
                    upstream_stop
                ),
                is_partial = alt_upstream_exon.donor < pas_overlap_junction.donor
            )
            bypass_events.insert(0, first_bypass_event)
        elif any(exon.donor <= upstream_stop < exon.acceptor for exon in (anchor_exons if is_deletion else other_exons)):
            downstream_exons[0] = evolve(downstream_exons[0], donor=upstream_stop+1)
            bypass_events[0] = evolve(bypass_events[0], exon=downstream_exons[0], is_partial=True)
        apa_event = APAEvent.from_basic_events(bypass_events)
    else:
        apa_event = None
    return splice_events, tss_event, apa_event

# %%
@frozen
class AlignmentBlock:
    category: 'TranscriptLevelAlignmentCategory' = field(init=False)
    anchor_range: range = field(factory=range, order=attrgetter('start', 'stop'))
    other_range: range = field(factory=range, order=attrgetter('start, stop'))
    
    def __attrs_post_init__(self):
        A, O = len(self.anchor_range), len(self.other_range)
        if A == O > 0:
            object.__setattr__(self, 'category', TranscriptLevelAlignmentCategory.MATCH)
        elif O > A == 0:
            object.__setattr__(self, 'category', TranscriptLevelAlignmentCategory.INSERTION)
        elif A > O == 0:
            object.__setattr__(self, 'category', TranscriptLevelAlignmentCategory.DELETION)
        else:
            raise ValueError(f'Invalid ranges {self.anchor_range} and {self.other_range}')
    
    def __repr__(self):
        return f'{self.category}({self.anchor_range.start}:{self.anchor_range.stop}|{self.other_range.start}:{self.other_range.stop})'


@define
class TranscriptAlignment:
    anchor: 'Transcript'
    other: 'Transcript' = field()
    events: tuple['CompoundTranscriptEvent'] = field(converter=sort_events)
    anchor_map: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    other_map: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    
    @other.validator
    def _check_transcripts(self, attribute, value):
        if value.gene_id != self.anchor.gene_id:
            raise ValueError(f'{self.anchor} and {value} are from different genes')
    @events.validator
    def _check_events(self, attribute, value):
        if sum(event.delta_nt for event in value) != self.other.length - self.anchor.length:
            raise ValueError(f'TranscriptEvent lengths do not add up correctly')
    
    def __attrs_post_init__(self):
        if any(i.data.is_insertion for i in self.anchor_map.all_intervals if isinstance(i.data, BasicTranscriptEvent)):
            raise ValueError
        if any(i.data.is_deletion for i in self.other_map.all_intervals if isinstance(i.data, BasicTranscriptEvent)):
            raise ValueError
        anchor_matches = {i.data for i in self.anchor_map.all_intervals if isinstance(i.data, AlignmentBlock) and i.data.category is TranscriptLevelAlignmentCategory.MATCH}
        other_matches = {i.data for i in self.other_map.all_intervals if isinstance(i.data, AlignmentBlock) and i.data.category is TranscriptLevelAlignmentCategory.MATCH}
        if anchor_matches != other_matches:
            raise ValueError(f'{anchor_matches} != {other_matches}')

    @classmethod
    def from_transcripts(cls, anchor: 'Transcript', other: 'Transcript'):
        splice_events, tss_event, apa_event = call_transcript_events(anchor, other)
        events = splice_events.copy()
        if tss_event:
            events.append(tss_event)
        if apa_event:
            events.append(apa_event)

        # map all events to transcript coordinates
        def get_transcript_interval(event: 'BasicTranscriptEvent'):
            transcript = anchor if event.is_deletion else other
            start = transcript.get_nucleotide_from_coordinate(event.start.coordinate).position
            stop = transcript.get_nucleotide_from_coordinate(event.stop.coordinate).position
            return Interval(start-1, stop, event)
        
        event_to_interval: dict['BasicTranscriptEvent', 'Interval'] = dict()
        basic_to_compound: dict['BasicTranscriptEvent', 'CompoundTranscriptEvent'] = dict()
        for compound_event in events:
            for event in compound_event.members:
                event_to_interval[event] = get_transcript_interval(event)
                basic_to_compound[event] = compound_event
        
        def get_compound_map(basic_map: 'IntervalTree'):
            compound_map = IntervalTree()
            for compound_event, intervals in groupby(sorted(basic_map.all_intervals), key=lambda i: basic_to_compound[i.data]):
                if len(compound_event.members) == 1:
                    continue
                compound_interval = IntervalTree(intervals)
                compound_interval.merge_neighbors()
                compound_map.update(i._replace(data=compound_event) for i in compound_interval.all_intervals)
            return compound_map

        insertions, deletions = partition(attrgetter('is_deletion'), event_to_interval.keys())
        anchor_basic = IntervalTree(event_to_interval[event] for event in deletions)
        other_basic = IntervalTree(event_to_interval[event] for event in insertions)
        anchor_compound = get_compound_map(anchor_basic)
        other_compound = get_compound_map(other_basic)

        # determine deletion and insertion block ranges
        del_ranges = anchor_basic.copy()
        del_ranges.merge_neighbors()
        ins_ranges = other_basic.copy()
        ins_ranges.merge_neighbors()

        # determine match block ranges
        del_blocks, ins_blocks, match_blocks = [], [], []

        def add_match_block(length, anchor_start, other_start):
            if length > 0:
                match_block = AlignmentBlock(
                    range(anchor_start, anchor_start + length),
                    range(other_start, other_start + length)
                )
                match_blocks.append(match_block)
                return anchor_start + length, other_start + length
            else:
                return anchor_start, other_start

        anchor_pos, other_pos = 0, 0
        sorted_del_ranges = deque(sorted(i[:2] for i in del_ranges))
        sorted_ins_ranges = deque(sorted(i[:2] for i in ins_ranges))
        while sorted_del_ranges or sorted_ins_ranges:
            to_next_del_block = (sorted_del_ranges[0][0] - anchor_pos) if sorted_del_ranges else float('inf')
            to_next_ins_block = (sorted_ins_ranges[0][0] - other_pos) if sorted_ins_ranges else float('inf')
            anchor_pos, other_pos = add_match_block(min(to_next_del_block, to_next_ins_block), anchor_pos, other_pos)
            if to_next_del_block <= to_next_ins_block:
                del_start, del_stop = sorted_del_ranges.popleft()
                del_blocks.append(AlignmentBlock(
                    range(del_start, del_stop),
                    range(other_pos, other_pos)
                ))
                anchor_pos = del_stop
            if to_next_ins_block <= to_next_del_block:
                ins_start, ins_stop = sorted_ins_ranges.popleft()
                ins_blocks.append(AlignmentBlock(
                    range(anchor_pos, anchor_pos),
                    range(ins_start, ins_stop)
                ))
                other_pos = ins_stop
        assert anchor.length - anchor_pos == other.length - other_pos, f'{anchor.length:=}, {anchor_pos:=}, {other.length:=}, {other_pos:=}'
        add_match_block(anchor.length - anchor_pos, anchor_pos, other_pos)
            
        anchor_blocks = IntervalTree.from_tuples((block.anchor_range.start, block.anchor_range.stop, block) for block in match_blocks + del_blocks)
        other_blocks = IntervalTree.from_tuples((block.other_range.start, block.other_range.stop, block) for block in match_blocks + ins_blocks)

        anchor_map = anchor_blocks.union(anchor_compound).union(anchor_basic)
        other_map = other_blocks.union(other_compound).union(other_basic)
        return cls(anchor, other, events, anchor_map, other_map)

# %%
t0 = time()
all_alns: dict[tuple[str, str], 'TranscriptAlignment'] = dict()
for a, o, _ in df.sort_values('anchor').itertuples(index=False):
    anchor: 'Transcript' = txs[a]
    other: 'Transcript' = txs[o]

    aln = TranscriptAlignment.from_transcripts(anchor, other)
    all_alns[a, o] = aln
t1 = time()

# %%
for (a, o), aln in all_alns.items():
    print(aln)
    for event in aln.events:
        print(f'\t{type(event).__name__}')
        for e in event.members:
            print(f'\t\t{e}')
    print(f'{a} events')
    for i in sorted(filter(lambda i: isinstance(i.data, TranscriptEvent), aln.anchor_map)):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t' + getattr(i.data, 'code', type(i.data).__name__))
    print(f'{a} blocks')
    for i in sorted(filter(lambda i: isinstance(i.data, AlignmentBlock), aln.anchor_map)):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t{i.data}')
    print(f'{o} events')
    for i in sorted(filter(lambda i: isinstance(i.data, TranscriptEvent), aln.other_map)):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t' + getattr(i.data, 'code', type(i.data).__name__))
    print(f'{o} blocks')
    for i in sorted(filter(lambda i: isinstance(i.data, AlignmentBlock), aln.other_map)):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t{i.data}')

# %%
for (a, o), aln in all_alns.items():
    event_sum = sum(event.delta_nt for event in aln.events)
    tx_diff = txs[o].length - txs[a].length
    if event_sum != tx_diff:
        print(f'{a}, {o}: {event_sum} != {tx_diff}')
        for event in aln.events:
            print(f'\t{type(event).__name__}: {event.delta_nt}')

# %%
EVENT_COLORS = {
    IntronSpliceEvent: '#e69138',
    DonorSpliceEvent: '#6aa84f',
    AcceptorSpliceEvent: '#674ea7',
    ExonSpliceEvent: '#3d85c6',
    ExonBypassEvent: '#bebebe',
}

for a, others in df.groupby('anchor')['other']:
    fig_path = f'../output/splice-test/{txs[a].gene.name}.png'
    if not os.path.isfile(fig_path):
        transcripts = [txs[a]] + [txs[other] for other in others]
        isoplot = IsoformPlot(transcripts)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        
        for i, other in enumerate(others, start=1):
            for event in all_alns[a, other].events:
                for base_event in event.members:
                    isoplot.draw_region(
                        i,
                        base_event.start.coordinate,
                        base_event.stop.coordinate,
                        y_offset = -0.5,
                        height = 0.1,
                        facecolor = EVENT_COLORS[type(base_event)]
                    )
        isoplot.fig.set_size_inches(10, 0.2 + 0.4 * len(transcripts))
        plt.savefig(fig_path, dpi=200, bbox_inches='tight', facecolor='w')
        print('saved '+fig_path)
        plt.close(isoplot.fig)

# %%
print(f'{t1 - t0:.3g} s total')
print(f'{(t1 - t0)/len(all_alns):.3g}s per alignment')

# %%
