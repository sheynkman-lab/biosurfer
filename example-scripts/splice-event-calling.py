# %%
import heapq
import os
import warnings
from abc import ABC, abstractmethod
from collections import deque
from collections.abc import Iterable
from functools import lru_cache
from itertools import chain, groupby
from operator import attrgetter, methodcaller
from time import time
from typing import Union
from warnings import warn

import numpy as np
import pandas as pd
from attrs import define, evolve, field, frozen, validators
from biosurfer.core.constants import AlignmentCategory
from biosurfer.core.constants import \
    ProteinLevelAlignmentCategory as SeqAlignCat
from biosurfer.core.constants import \
    TranscriptLevelAlignmentCategory as CodonAlignCat
from biosurfer.core.database import Database
from biosurfer.core.helpers import (BisectDict, ExceptionLogger, Interval,
                                    IntervalTree, get_interval_overlap_graph)
from biosurfer.core.models.biomolecules import Exon, Transcript, Protein
from biosurfer.core.models.nonpersistent import Codon, Junction, Position
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
            raise ValueError(f'{value}')
        if len({event.is_deletion for event in value[:-1]}) > 1:
            raise ValueError(f'{value}')


@frozen(eq=True)
class APAEvent(CompoundTranscriptEvent):
    members: tuple['ExonBypassEvent', ...] = field(converter=sort_events)

    @members.validator
    def _check_members(self, attribute, value: tuple['ExonBypassEvent', ...]):
        if any(event.is_partial for event in value[1:]):
            raise ValueError(f'{value}')
        if len({event.is_deletion for event in value[1:]}) > 1:
            raise ValueError(f'{value}')


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
CACHE_SIZE = 2**8


@frozen(order=True)
class AlignmentBlock(ABC):
    anchor_range: range = field(factory=range, order=attrgetter('start', 'stop'))
    other_range: range = field(factory=range, order=attrgetter('start', 'stop'))
    category: 'AlignmentCategory' = field(default=None, order=False)
    
    @abstractmethod
    def __attrs_post_init__(self):
        pass
    
    def __repr__(self):
        return f'{self.category}({self.anchor_range.start}:{self.anchor_range.stop}|{self.other_range.start}:{self.other_range.stop})'


@frozen(order=True, repr=False)
class TranscriptAlignmentBlock(AlignmentBlock):
    category: 'SeqAlignCat' = field(init=False)

    def __attrs_post_init__(self):
        A, O = len(self.anchor_range), len(self.other_range)
        if A == O > 0:
            object.__setattr__(self, 'category', SeqAlignCat.MATCH)
        elif O > A == 0:
            object.__setattr__(self, 'category', SeqAlignCat.INSERTION)
        elif A > O == 0:
            object.__setattr__(self, 'category', SeqAlignCat.DELETION)
        else:
            raise ValueError(f'Invalid ranges {self.anchor_range} and {self.other_range}')


@frozen(repr=False)
class CodonAlignmentBlock(AlignmentBlock):
    category: 'CodonAlignCat' = field(kw_only=True)

    def __attrs_post_init__(self):
        A, O = len(self.anchor_range), len(self.other_range)
        if A == O == 0:
            raise ValueError(f'Invalid ranges {self.anchor_range} and {self.other_range}')


@define
class TranscriptAlignment:
    anchor: 'Transcript'
    other: 'Transcript' = field()
    events: tuple['CompoundTranscriptEvent'] = field(converter=sort_events)
    anchor_events: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    anchor_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    other_events: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    other_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False)

    @other.validator
    def _check_transcripts(self, attribute, value):
        if value.gene_id != self.anchor.gene_id:
            raise ValueError(f'{self.anchor} and {value} are from different genes')
    
    def __attrs_post_init__(self):
        if any(i.data.is_insertion for i in self.anchor_events if isinstance(i.data, BasicTranscriptEvent)):
            raise ValueError
        if any(i.data.is_deletion for i in self.other_events if isinstance(i.data, BasicTranscriptEvent)):
            raise ValueError
        anchor_matches = {i.data for i in self.anchor_blocks if i.data.category is SeqAlignCat.MATCH}
        other_matches = {i.data for i in self.other_blocks if i.data.category is SeqAlignCat.MATCH}
        if anchor_matches != other_matches:
            raise ValueError(f'{anchor_matches} != {other_matches}')
        total_delta_nt = sum(event.delta_nt for event in self.events)
        tx_length_diff = self.other.length - self.anchor.length
        if total_delta_nt != tx_length_diff:
            raise ValueError(f'TranscriptEvent lengths add up to {total_delta_nt}; expected {tx_length_diff}')

    @property
    def basic_events(self) -> tuple['BasicTranscriptEvent', ...]:
        return tuple(chain.from_iterable(event.members for event in self.events))
    
    @property
    def blocks(self) -> tuple['TranscriptAlignmentBlock', ...]:
        return tuple(sorted({i.data for i in chain(self.anchor_blocks, self.other_blocks)}))

    @classmethod
    @lru_cache(maxsize=CACHE_SIZE)
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
        anchor_events = anchor_compound.union(anchor_basic)
        other_events = other_compound.union(other_basic)

        # determine deletion and insertion block ranges
        del_ranges = anchor_basic.copy()
        del_ranges.merge_neighbors()
        ins_ranges = other_basic.copy()
        ins_ranges.merge_neighbors()

        # determine match block ranges
        blocks = []
        position = {'anchor': 0, 'other': 0}
        sorted_del_ranges = deque(sorted(i[:2] for i in del_ranges))
        sorted_ins_ranges = deque(sorted(i[:2] for i in ins_ranges))

        def add_match_block(length: int):
            if length > 0:
                match_block = TranscriptAlignmentBlock(
                    range(position['anchor'], position['anchor'] + length),
                    range(position['other'], position['other'] + length)
                )
                blocks.append(match_block)
            position['anchor'] += length
            position['other'] += length

        def add_del_block():
            del_start, del_stop = sorted_del_ranges.popleft()
            blocks.append(TranscriptAlignmentBlock(
                range(del_start, del_stop),
                range(position['other'], position['other'])
            ))
            position['anchor'] = del_stop
        
        def add_ins_block():
            ins_start, ins_stop = sorted_ins_ranges.popleft()
            blocks.append(TranscriptAlignmentBlock(
                range(position['anchor'], position['anchor']),
                range(ins_start, ins_stop)
            ))
            position['other'] = ins_stop

        while sorted_del_ranges or sorted_ins_ranges:
            to_next_del_block = (sorted_del_ranges[0][0] - position['anchor']) if sorted_del_ranges else float('inf')
            to_next_ins_block = (sorted_ins_ranges[0][0] - position['other']) if sorted_ins_ranges else float('inf')
            add_match_block(min(to_next_del_block, to_next_ins_block))
            if to_next_del_block < to_next_ins_block:
                add_del_block()
            elif to_next_ins_block < to_next_del_block:
                add_ins_block()
            else:
                anchor_pos_genomic = anchor.get_genome_coord_from_transcript_coord(position['anchor'])
                other_pos_genomic = other.get_genome_coord_from_transcript_coord(position['other'])
                if anchor_pos_genomic < other_pos_genomic:
                    add_del_block()
                    add_ins_block()
                else:
                    add_ins_block()
                    add_del_block()
        assert anchor.length - position['anchor'] == other.length - position['other'], f'{position=}, {anchor.length=}, {other.length=}'
        add_match_block(anchor.length - position['anchor'])
            
        anchor_blocks = IntervalTree.from_tuples((block.anchor_range.start, block.anchor_range.stop, block) for block in blocks if block.anchor_range)
        other_blocks = IntervalTree.from_tuples((block.other_range.start, block.other_range.stop, block) for block in blocks if block.other_range)

        return cls(anchor, other, events, anchor_events, anchor_blocks, other_events, other_blocks)


@define
class ProteinAlignment:
    anchor: 'Protein'
    other: 'Protein'
    anchor_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    other_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False)

    @classmethod
    @lru_cache(maxsize=CACHE_SIZE)
    def from_proteins(cls, anchor: 'Protein', other: 'Protein'):
        tx_aln = TranscriptAlignment.from_transcripts(anchor.transcript, other.transcript)
        anchor_orf_range = range(anchor.orf.transcript_start - 1, anchor.orf.transcript_stop)
        other_orf_range = range(other.orf.transcript_start - 1, other.orf.transcript_stop)
        anchor_orf_len = len(anchor_orf_range)
        other_orf_len = len(other_orf_range)
        anchor_pr_len = anchor.length
        other_pr_len = other.length

        def compare_ranges(a: range, b: range):
            if a.start < b.stop and b.start < a.stop:
                return 0
            elif a.stop <= b.start:
                return -1
            else:
                return 1

        # convert transcript-relative coords to ORF-relative coords
        tx_blocks = deque(
            (
                block.category,
                range(block.anchor_range.start - anchor_orf_range.start, block.anchor_range.stop - anchor_orf_range.start),
                range(block.other_range.start - other_orf_range.start, block.other_range.stop - other_orf_range.start)
            ) for block in tx_aln.blocks
        )

        # determine whether N-termini are different
        # anchor_ic_coord = anchor.transcript.get_genome_coord_from_transcript_coord(anchor_orf_range[0])
        # other_ic_coord = other.transcript.get_genome_coord_from_transcript_coord(other_orf_range[0])
        # if anchor_ic_coord != other_ic_coord:
        #     pass

        anchor_cd_blocks, other_cd_blocks = [], []
        frame = 0
        frame_to_category = {
            0: CodonAlignCat.MATCH,
            1: CodonAlignCat.FRAME_AHEAD,
            2: CodonAlignCat.FRAME_BEHIND
        }
        while tx_blocks:
            tx_category, anchor_tx_range, other_tx_range = tx_blocks.popleft()

            # if block overlaps an ORF boundary, split it up
            def split_paired_ranges(a: range, b: range, a_split: int):
                a_left, a_right = range(a.start, a_split), range(a_split, a.stop)
                b_split = b.start + len(a_left)
                b_left, b_right = range(b.start, b_split), range(b_split, b.stop)
                return a_left, a_right, b_left, b_right
            if anchor_tx_range.start < 0 < anchor_tx_range.stop:
                anchor_tx_range, next_anchor_range, other_tx_range, next_other_range = split_paired_ranges(anchor_tx_range, other_tx_range, 0)
                tx_blocks.appendleft((tx_category, next_anchor_range, next_other_range))
            if other_tx_range.start < 0 < other_tx_range.stop:
                other_tx_range, next_other_range, anchor_tx_range, next_anchor_range = split_paired_ranges(other_tx_range, anchor_tx_range, 0)
                tx_blocks.appendleft((tx_category, next_anchor_range, next_other_range))
            if anchor_tx_range.start < anchor_orf_len < anchor_tx_range.stop:
                anchor_tx_range, next_anchor_range, other_tx_range, next_other_range = split_paired_ranges(anchor_tx_range, other_tx_range, anchor_orf_len)
                tx_blocks.appendleft((tx_category, next_anchor_range, next_other_range))
            if other_tx_range.start < other_orf_len < other_tx_range.stop:
                other_tx_range, next_other_range, anchor_tx_range, next_anchor_range = split_paired_ranges(other_tx_range, anchor_tx_range, other_orf_len)
                tx_blocks.appendleft((tx_category, next_anchor_range, next_other_range))
            
            
            # skip blocks that are outside both ORF ranges
            outside_anchor_orf = compare_ranges(anchor_tx_range, range(0, len(anchor_orf_range)))
            outside_other_orf = compare_ranges(other_tx_range, range(0, len(other_orf_range)))
            if outside_anchor_orf and outside_other_orf:
                continue
            
            # convert block range to protein coords
            if outside_anchor_orf < 0:
                anchor_pr_range = range(0, 0)
            elif outside_anchor_orf > 0:
                anchor_pr_range = range(anchor_pr_len, anchor_pr_len)
            else:
                anchor_pr_range = range(anchor_tx_range.start//3, anchor_tx_range.stop//3)
            if outside_other_orf < 0:
                other_pr_range = range(0, 0)
            elif outside_other_orf > 0:
                other_pr_range = range(other_pr_len, other_pr_len)
            else:
                other_pr_range = range(other_tx_range.start//3, other_tx_range.stop//3)

            # infer codon block categories
            if tx_category is SeqAlignCat.MATCH:
                if not anchor_pr_range:
                    cd_category = CodonAlignCat.TRANSLATED
                elif not other_pr_range:
                    cd_category = CodonAlignCat.UNTRANSLATED
                else:
                    cd_category = frame_to_category[frame]
            else:
                if tx_category is SeqAlignCat.DELETION:
                    frame = (frame + len(anchor_tx_range)) % 3
                    cd_category = CodonAlignCat.DELETION
                elif tx_category is SeqAlignCat.INSERTION:
                    frame = (frame + len(other_tx_range)) % 3
                    cd_category = CodonAlignCat.INSERTION
                else:
                    raise RuntimeError
            if anchor_pr_range or other_pr_range:
                cd_block = CodonAlignmentBlock(
                    anchor_pr_range,
                    other_pr_range,
                    category = cd_category
                )
                if anchor_pr_range:
                    anchor_cd_blocks.append(cd_block)
                if other_pr_range:
                    other_cd_blocks.append(cd_block)
        
        anchor_blocks = IntervalTree.from_tuples(
            (block.anchor_range.start, block.anchor_range.stop, block)
            for block in anchor_cd_blocks
        )
        other_blocks = IntervalTree.from_tuples(
            (block.other_range.start, block.other_range.stop, block)
            for block in other_cd_blocks
        )
        return cls(anchor, other, anchor_blocks, other_blocks)


# %%
all_alns: dict[tuple[str, str], tuple['TranscriptAlignment', 'ProteinAlignment']] = dict()
t0 = time()
for a, o, _ in df.sort_values('anchor').itertuples(index=False):
    anchor: 'Transcript' = txs[a]
    other: 'Transcript' = txs[o]

    pr_aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
    tx_aln = TranscriptAlignment.from_transcripts(anchor, other)
    all_alns[a, o] = tx_aln, pr_aln
t1 = time()

# %%
for (a, o), (tx_aln, pr_aln) in all_alns.items():
    print(tx_aln)
    for event in tx_aln.events:
        print(f'\t{type(event).__name__}')
        for e in event.members:
            print(f'\t\t{e}')
    print(f'{a} events')
    for i in sorted(tx_aln.anchor_events):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t' + getattr(i.data, 'code', type(i.data).__name__))
    print(f'{a} tx blocks')
    for i in sorted(tx_aln.anchor_blocks):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t{i.data}')
    print(f'{a} pr blocks')
    for i in sorted(pr_aln.anchor_blocks):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t{i.data}')
    print(f'{o} events')
    for i in sorted(tx_aln.other_events):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t' + getattr(i.data, 'code', type(i.data).__name__))
    print(f'{o} tx blocks')
    for i in sorted(tx_aln.other_blocks):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t{i.data}')
    print(f'{o} pr blocks')
    for i in sorted(pr_aln.other_blocks):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t{i.data}')

# %%
# EVENT_COLORS = {
#     IntronSpliceEvent: '#e69138',
#     DonorSpliceEvent: '#6aa84f',
#     AcceptorSpliceEvent: '#674ea7',
#     ExonSpliceEvent: '#3d85c6',
#     ExonBypassEvent: '#bebebe',
# }

# for a, others in df.groupby('anchor')['other']:
#     fig_path = f'../output/splice-test/{txs[a].gene.name}.png'
#     if not os.path.isfile(fig_path):
#         transcripts = [txs[a]] + [txs[other] for other in others]
#         isoplot = IsoformPlot(transcripts)
#         isoplot.draw_all_isoforms()
#         isoplot.draw_frameshifts()
        
#         for i, other in enumerate(others, start=1):
#             for event in all_alns[a, other].events:
#                 for base_event in event.members:
#                     isoplot.draw_region(
#                         i,
#                         base_event.start.coordinate,
#                         base_event.stop.coordinate,
#                         y_offset = -0.5,
#                         height = 0.1,
#                         facecolor = EVENT_COLORS[type(base_event)]
#                     )
#         isoplot.fig.set_size_inches(10, 0.2 + 0.4 * len(transcripts))
#         plt.savefig(fig_path, dpi=200, bbox_inches='tight', facecolor='w')
#         print('saved '+fig_path)
#         plt.close(isoplot.fig)

# %%
print(f'{t1 - t0:.3g} s total')
print(f'{(t1 - t0)/len(all_alns):.3g}s per alignment')

# %%
