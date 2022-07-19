import heapq
import warnings
from abc import ABC, abstractmethod
from operator import attrgetter, methodcaller
from typing import Iterable, Union

from attrs import evolve, field, frozen
from biosurfer.core.helpers import get_interval_overlap_graph
from biosurfer.core.models.biomolecules import Exon, Transcript
from biosurfer.core.models.nonpersistent import GenomeRange, Position, Junction
from graph_tool import GraphView
from graph_tool.topology import is_bipartite, shortest_path, label_components
from more_itertools import first, windowed


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

@frozen(eq=True)
class ExonBypassEvent(BasicTranscriptEvent):
    exon: 'GenomeRange'  # TODO: replace with updated Exon object
    is_partial: bool = False

    @property
    def start(self) -> 'Position':
        return self.exon.begin
    
    @property
    def stop(self) -> 'Position':
        return self.exon.end


EVENT_CODES = {
    IntronSpliceEvent: ('I', 'i'),
    DonorSpliceEvent: ('D', 'd'),
    AcceptorSpliceEvent: ('A', 'a'),
    ExonSpliceEvent: ('E', 'e'),
    ExonBypassEvent: ('B', 'b')
}


def get_event_code(events: Iterable['BasicTranscriptEvent']):
    # return ''.join(EVENT_CODES[type(event)][event.is_deletion] for event in events)
    code = ''
    for event in events:
        if getattr(event, 'is_partial', False):
            code += 'p' if event.is_deletion else 'P'
        else:
            code += EVENT_CODES[type(event)][event.is_deletion]
    return code


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
        code = get_event_code(events)
        anchor_junctions = sorted({junc for event in events for junc in event.anchor_junctions}, key=attrgetter('donor'))
        other_junctions = sorted({junc for event in events for junc in event.other_junctions}, key=attrgetter('donor'))
        return cls(members=events, code=code, anchor_junctions=tuple(anchor_junctions), other_junctions=tuple(other_junctions))


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
    chr = anchor.gene.chromosome_id
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
        warnings.warn(f'Overlap graph not bipartite for {anchor} | {other}')
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
        return GenomeRange(*sorted((Position(chr, strand, exon_obj.start), Position(chr, strand, exon_obj.stop))))
    
    anchor_exons = [get_exon(exon) for exon in anchor.exons]
    other_exons = [get_exon(exon) for exon in other.exons]
    upstream_exons = sorted((exon for exon in set(anchor_exons) | set(other_exons) if exon.begin < downstream_start), key=attrgetter('begin'))
    downstream_exons = sorted((exon for exon in set(anchor_exons) | set(other_exons) if upstream_stop < exon.end), key=attrgetter('begin'))
    # call TSS event
    if upstream_exons:
        is_deletion = downstream_start == other_start
        bypass_events = [ExonBypassEvent(is_deletion, exon) for exon in upstream_exons]
        if tss_overlap_junction:
            alt_downstream_exon = (other_exons if is_deletion else anchor_exons)[0]
            last_bypass_event = ExonBypassEvent(
                not is_deletion,
                GenomeRange(
                    downstream_start,
                    min(alt_downstream_exon.end, tss_overlap_junction.acceptor)
                ),
                is_partial = tss_overlap_junction.acceptor < alt_downstream_exon.end
            )
            bypass_events.append(last_bypass_event)
        elif any(exon.begin < downstream_start <= exon.end for exon in (anchor_exons if is_deletion else other_exons)):
            upstream_exons[-1] = evolve(upstream_exons[-1], end=downstream_start-1)
            bypass_events[-1] = evolve(bypass_events[-1], exon=upstream_exons[-1], is_partial=True)
        tss_event = TSSEvent.from_basic_events(bypass_events)
    else:
        tss_event = None
    # call APA event
    if downstream_exons:
        is_deletion = upstream_stop == other_stop
        bypass_events = [ExonBypassEvent(is_deletion, exon) for exon in downstream_exons]
        if pas_overlap_junction:
            alt_upstream_exon = (other_exons if is_deletion else anchor_exons)[-1]
            first_bypass_event = ExonBypassEvent(
                not is_deletion,
                GenomeRange(
                    max(alt_upstream_exon.begin, pas_overlap_junction.donor),
                    upstream_stop
                ),
                is_partial = alt_upstream_exon.begin < pas_overlap_junction.donor
            )
            bypass_events.insert(0, first_bypass_event)
        elif any(exon.begin <= upstream_stop < exon.end for exon in (anchor_exons if is_deletion else other_exons)):
            downstream_exons[0] = evolve(downstream_exons[0], begin=upstream_stop+1)
            bypass_events[0] = evolve(bypass_events[0], exon=downstream_exons[0], is_partial=True)
        apa_event = APAEvent.from_basic_events(bypass_events)
    else:
        apa_event = None
    return splice_events, tss_event, apa_event
