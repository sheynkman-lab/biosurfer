from abc import ABC, abstractmethod
from collections import deque
from functools import lru_cache
from itertools import chain, groupby, tee
from operator import attrgetter, itemgetter
from typing import TYPE_CHECKING

from attrs import define, field, frozen
from biosurfer.core.constants import \
    ProteinLevelAlignmentCategory as SeqAlignCat
from biosurfer.core.constants import \
    TranscriptLevelAlignmentCategory as CodonAlignCat
from biosurfer.core.helpers import Interval, IntervalTree
from biosurfer.core.models.biomolecules import Protein, Transcript
from biosurfer.core.splice_events import (BasicTranscriptEvent,
                                          call_transcript_events, sort_events)
from more_itertools import first, last, partition

if TYPE_CHECKING:
    from biosurfer.core.constants import AlignmentCategory
    from biosurfer.core.splice_events import CompoundTranscriptEvent


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
    events: tuple['CompoundTranscriptEvent', ...] = field(converter=sort_events)
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

    @property
    def blocks(self) -> tuple['CodonAlignmentBlock', ...]:
        return tuple(sorted({i.data for i in chain(self.anchor_blocks, self.other_blocks)}))

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

        def split_paired_ranges(a: range, b: range, a_split: int):
            assert a_split in a
            a_left, a_right = range(a.start, a_split), range(a_split, a.stop)
            if b:
                b_split = b.start + len(a_left)
                b_left, b_right = range(b.start, b_split), range(b_split, b.stop)
            else:
                b_left, b_right = b, b
            return a_left, a_right, b_left, b_right
        
        cd_block_ranges = []
        frame_to_category = {
            0: CodonAlignCat.MATCH,
            1: CodonAlignCat.FRAME_AHEAD,
            2: CodonAlignCat.FRAME_BEHIND
        }
        while tx_blocks:
            tx_category, anchor_tx_range, other_tx_range = tx_blocks.popleft()

            # if block overlaps an ORF boundary, split it up
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
            anchor_pr_start = anchor_tx_range.start//3
            anchor_pr_stop = anchor_tx_range.stop//3
            other_pr_start = other_tx_range.start//3
            other_pr_stop = other_tx_range.stop//3
            if outside_anchor_orf < 0:
                anchor_pr_range = range(0, 0)
            elif outside_anchor_orf > 0:
                anchor_pr_range = range(anchor_pr_len, anchor_pr_len)
            else:
                anchor_pr_start = anchor_tx_range.start//3
                anchor_pr_range = range(anchor_pr_start, anchor_pr_stop)
            if outside_other_orf < 0:
                other_pr_range = range(0, 0)
            elif outside_other_orf > 0:
                other_pr_range = range(other_pr_len, other_pr_len)
            else:
                other_pr_start = other_tx_range.start//3
                other_pr_range = range(other_pr_start, other_pr_stop)

            if not (anchor_pr_range or other_pr_range):
                continue

            # infer categories of codon blocks
            frameshift = (other_tx_range.start - anchor_tx_range.start) % 3
            edge_first, edge_last = False, False
            if tx_category is SeqAlignCat.MATCH:
                if not anchor_pr_range:
                    cd_category = CodonAlignCat.TRANSLATED
                elif not other_pr_range:
                    cd_category = CodonAlignCat.UNTRANSLATED
                else:
                    edge_first = frameshift == 0 and anchor_tx_range.start % 3 == 1
                    edge_last = frameshift == 0 and anchor_tx_range.stop % 3 == 2
                    cd_category = frame_to_category[frameshift]
                    # TODO: detect triple overlap when len(anchor_pr_range) != len(other_pr_range)
            elif tx_category is SeqAlignCat.DELETION:
                cd_category = CodonAlignCat.DELETION
            elif tx_category is SeqAlignCat.INSERTION:
                cd_category = CodonAlignCat.INSERTION
            else:
                raise RuntimeError
            if edge_first:
                anchor_edge_first, anchor_pr_range = anchor_pr_range[:1], anchor_pr_range[1:]
                other_edge_first, other_pr_range = other_pr_range[:1], other_pr_range[1:]
                edge_first_category = CodonAlignCat.EDGE_MATCH if anchor.residues[anchor_edge_first[0]].amino_acid == other.residues[other_edge_first[0]].amino_acid else CodonAlignCat.EDGE_MISMATCH
            if edge_last:
                anchor_pr_range, anchor_edge_last = anchor_pr_range[:-1], anchor_pr_range[-1:]
                other_pr_range, other_edge_last = other_pr_range[:-1], other_pr_range[-1:]
                edge_last_category = CodonAlignCat.EDGE_MATCH if anchor.residues[anchor_edge_last[0]].amino_acid == other.residues[other_edge_last[0]].amino_acid else CodonAlignCat.EDGE_MISMATCH
            if edge_first:
                cd_block_ranges.append((edge_first_category, anchor_edge_first, other_edge_first))
            cd_block_ranges.append((cd_category, anchor_pr_range, other_pr_range))
            if edge_last:
                cd_block_ranges.append((edge_last_category, anchor_edge_last, other_edge_last))
        
        # merge consecutive codon blocks w/ same category
        cd_blocks = []
        for category, block_group in groupby(cd_block_ranges, key=itemgetter(0)):
            if category is CodonAlignCat.MATCH:
                cd_blocks.extend(CodonAlignmentBlock(block[1], block[2], category=category) for block in block_group)
            else:
                g1, g2 = tee(block_group)
                first_block = first(g1)
                last_block = last(g2)
                anchor_start = first_block[1].start
                anchor_stop = last_block[1].stop
                other_start = first_block[2].start
                other_stop = last_block[2].stop
                cd_blocks.append(CodonAlignmentBlock(range(anchor_start, anchor_stop), range(other_start, other_stop), category=category))
        
        for block in cd_blocks:
            if block.category not in {CodonAlignCat.EDGE_MATCH, CodonAlignCat.EDGE_MISMATCH}:
                continue
            anchor_edge_res = anchor.residues[block.anchor_range[0]]
            other_edge_res = other.residues[block.other_range[0]]
            assert anchor_edge_res.codon[2:] == other_edge_res.codon[2:] or anchor_edge_res.codon[:2] == other_edge_res.codon[:2], f'{anchor_edge_res}, {anchor_edge_res.codon} | {other_edge_res}, {other_edge_res.codon}'

        anchor_blocks = IntervalTree.from_tuples(
            (block.anchor_range.start, block.anchor_range.stop, block)
            for block in cd_blocks if block.anchor_range
        )
        other_blocks = IntervalTree.from_tuples(
            (block.other_range.start, block.other_range.stop, block)
            for block in cd_blocks if block.other_range
        )
        return cls(anchor, other, anchor_blocks, other_blocks)
