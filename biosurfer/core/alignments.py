from abc import ABC, abstractmethod
from collections import deque
from functools import cached_property, lru_cache
from itertools import chain, groupby, tee
from operator import attrgetter, itemgetter
from typing import TYPE_CHECKING
import os
from attrs import define, evolve, field, frozen
from biosurfer.core.constants import ANCHOR_EXCLUSIVE, FRAMESHIFT, CD_DEL_INS, OTHER_EXCLUSIVE, SEQ_DEL_INS, SPLIT_CODON, CodonAlignmentCategory as CodonAlignCat
from biosurfer.core.constants import SequenceAlignmentCategory as SeqAlignCat
from biosurfer.core.helpers import Interval, IntervalTree
from biosurfer.core.models.biomolecules import Protein, Transcript
from biosurfer.core.models.features import ProjectedFeature, ProteinFeature
from biosurfer.core.splice_events import (BasicTranscriptEvent, TranscriptEvent,
                                          call_transcript_events, sort_events)
from more_itertools import first, last, one, only, partition, windowed

if TYPE_CHECKING:
    from biosurfer.core.constants import AlignmentCategory
    from biosurfer.core.splice_events import CompoundTranscriptEvent


CACHE_SIZE = 2**8

def check_block_ranges(instance, attribute, value: 'IntervalTree'):
    starts = sorted(i.begin for i in value)
    stops = sorted(i.end for i in value)
    if starts[0] < 0:
        raise ValueError(f'Block ranges cannot be negative')
    if starts[0] > 0:
        raise ValueError(f'Block ranges do not cover (0, {starts[0]})')
    for i, (start, stop) in enumerate(zip(starts[1:], stops[:-1])):
        if start != stop:
            raise ValueError(f'Gap or overlap between block ranges ({starts[i]}, {stop}) and ({start}, {stops[i]})')


@frozen(order=True)
class AlignmentBlock(ABC):
    anchor_range: range = field(factory=range, order=attrgetter('start', 'stop'))
    other_range: range = field(factory=range, order=attrgetter('start', 'stop'))
    category: 'AlignmentCategory' = field(default=None, order=False)
    
    def __attrs_post_init__(self):
        A, O = len(self.anchor_range), len(self.other_range)
        if A == O == 0:
            raise ValueError(f'Invalid ranges {self.anchor_range} and {self.other_range}')

    def __repr__(self):
        return f'{self.category}({self.anchor_range.start}:{self.anchor_range.stop}|{self.other_range.start}:{self.other_range.stop})'
    
    @property
    def delta_length(self):
        return len(self.other_range) - len(self.anchor_range)

    def project_coordinate(self, coord: int, *, from_anchor: bool = True):
        source, target = (self.anchor_range, self.other_range) if from_anchor else (self.other_range, self.anchor_range)
        index = source.index(coord)
        try:
            return target[index]
        except IndexError:
            return None


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


@frozen(order=True, repr=False)
class ProteinAlignmentBlock(AlignmentBlock):
    category: 'SeqAlignCat' = field(kw_only=True)
    ragged5: bool = field(default=False, order=False)
    ragged3: bool = field(default=False, order=False)

    def __attrs_post_init__(self):
        super().__attrs_post_init__()
        if self.ragged and self.category is SeqAlignCat.MATCH:
            raise ValueError(f'Match protein alignment blocks cannot be ragged')
    
    @property
    def ragged(self):
        return self.ragged5 or self.ragged3


class ProjectionMixin:
    def range_to_blocks(self, start: int, stop: int, *, from_anchor: bool = True):
        mapper: 'IntervalTree' = self.anchor_blocks if from_anchor else self.other_blocks
        blocks: list['AlignmentBlock'] = [interval.data for interval in sorted(mapper.overlap(start, stop))]
        if not blocks:
            raise ValueError(f'Could not locate mapping in {self.anchor if from_anchor else self.other} for range({start}, {stop})')
        first_block = blocks[0]
        last_block = blocks[-1]
        first_block_source = first_block.anchor_range if from_anchor else first_block.other_range
        last_block_source = last_block.anchor_range if from_anchor else last_block.other_range
        first_block_offset = start - first_block_source.start
        last_block_offset = stop - last_block_source.start
        if first_block is last_block:
            blocks[0] = evolve(
                first_block,
                anchor_range = first_block.anchor_range[first_block_offset:last_block_offset],
                other_range = first_block.other_range[first_block_offset:last_block_offset]
            )
        else:
            blocks[0] = evolve(
                first_block,
                anchor_range = first_block.anchor_range[first_block_offset:],
                other_range = first_block.other_range[first_block_offset:]
            )
            blocks[-1] = evolve(
                last_block,
                anchor_range = last_block.anchor_range[:last_block_offset],
                other_range = last_block.other_range[:last_block_offset]
            )
        return blocks

    def project_range(self, start: int, stop: int, *, from_anchor: bool = True) -> range:
        blocks = self.range_to_blocks(start, stop, from_anchor=from_anchor)
        if from_anchor:
            return range(blocks[0].other_range.start, blocks[-1].other_range.stop)
        else:
            return range(blocks[0].anchor_range.start, blocks[-1].anchor_range.stop)
    
    def project_coordinate(self, coord: int, *, from_anchor: bool = True):
        mapper: 'IntervalTree' = self.anchor_blocks if from_anchor else self.other_blocks
        block: 'AlignmentBlock' = one(mapper.at(coord)).data
        return block.project_coordinate(coord, from_anchor=from_anchor)


@define(eq=False)
class TranscriptAlignment(ProjectionMixin):
    anchor: 'Transcript'
    other: 'Transcript' = field()
    events: tuple['CompoundTranscriptEvent', ...] = field(converter=sort_events, repr=False)
    anchor_events: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    anchor_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False, validator=check_block_ranges)
    other_events: 'IntervalTree' = field(factory=IntervalTree, repr=False)
    other_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False, validator=check_block_ranges)
    event_to_block: dict['BasicTranscriptEvent', 'TranscriptAlignmentBlock'] = field(factory=dict, repr=False)
    block_to_events: dict['TranscriptAlignmentBlock', tuple['BasicTranscriptEvent', ...]] = field(factory=dict, repr=False)

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
        for event in self.basic_events:
            block = self.event_to_block.get(event, None)
            if event not in self.block_to_events.get(block, ()):
                raise ValueError(f'Broken event-to-block mapping for {event}')

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
            start = transcript.get_transcript_coord_from_genome_coord(event.start)
            stop = transcript.get_transcript_coord_from_genome_coord(event.stop) + 1
            return Interval(start, stop, event)
        
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
        del_ranges.merge_neighbors(data_reducer=lambda a, b: a + (b,), data_initializer=())
        ins_ranges = other_basic.copy()
        ins_ranges.merge_neighbors(data_reducer=lambda a, b: a + (b,), data_initializer=())

        # determine match block ranges
        blocks = []
        block_to_events = dict()
        event_to_block = dict()
        position = {'anchor': 0, 'other': 0}
        sorted_del_ranges = deque(sorted(map(tuple, del_ranges)))
        sorted_ins_ranges = deque(sorted(map(tuple, ins_ranges)))

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
            del_start, del_stop, events = sorted_del_ranges.popleft()
            del_block = TranscriptAlignmentBlock(
                range(del_start, del_stop),
                range(position['other'], position['other'])
            )
            blocks.append(del_block)
            block_to_events[del_block] = events
            for event in events:
                event_to_block[event] = del_block
            position['anchor'] = del_stop
        
        def add_ins_block():
            ins_start, ins_stop, events = sorted_ins_ranges.popleft()
            ins_block = TranscriptAlignmentBlock(
                range(position['anchor'], position['anchor']),
                range(ins_start, ins_stop)
            )
            blocks.append(ins_block)
            block_to_events[ins_block] = events
            for event in events:
                event_to_block[event] = ins_block
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
        assert anchor.length - position['anchor'] == other.length - position['other'], f'{position=}, {anchor.length=}, {other.length=}, {anchor=}'
        add_match_block(anchor.length - position['anchor'])
            
        anchor_blocks = IntervalTree.from_tuples((block.anchor_range.start, block.anchor_range.stop, block) for block in blocks if block.anchor_range)
        other_blocks = IntervalTree.from_tuples((block.other_range.start, block.other_range.stop, block) for block in blocks if block.other_range)

        return cls(anchor, other, events, anchor_events, anchor_blocks, other_events, other_blocks, event_to_block, block_to_events)


@define(eq=False)
class CodonAlignment(ProjectionMixin):
    anchor: 'Protein'
    other: 'Protein'
    anchor_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False, validator=check_block_ranges)
    other_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False, validator=check_block_ranges)
    tblock_to_cblocks: dict['TranscriptAlignmentBlock', tuple['CodonAlignmentBlock', ...]] = field(factory=dict, repr=False)
    cblock_to_tblock: dict['CodonAlignmentBlock', 'TranscriptAlignmentBlock'] = field(factory=dict, repr=False)

    def __attrs_post_init__(self):
        anchor_shared = {i.data for i in self.anchor_blocks if i.data.category not in ANCHOR_EXCLUSIVE}
        other_shared = {i.data for i in self.other_blocks if i.data.category not in OTHER_EXCLUSIVE}
        if anchor_shared != other_shared:
            raise ValueError(f'{anchor_shared} != {other_shared}')
        for cblock in filter(lambda cblock: cblock.category is not CodonAlignCat.MATCH, self.blocks):
            tblock = self.cblock_to_tblock.get(cblock, None)
            if tblock and cblock not in self.tblock_to_cblocks.get(tblock, ()):
                raise ValueError(f'Broken tblock-to-cblock mapping for {cblock}')

    @property
    def blocks(self) -> tuple['CodonAlignmentBlock', ...]:
        return tuple(sorted({i.data for i in chain(self.anchor_blocks, self.other_blocks)}))

    def project_feature(self, anchor_feature: 'ProteinFeature'):
        if anchor_feature.protein is not self.anchor:
            raise ValueError(f'{anchor_feature} is not a feature of {self.anchor}')
        other_range = self.project_range(anchor_feature.protein_start - 1, anchor_feature.protein_stop)
        if not other_range:
            return None
        other_blocks = self.range_to_blocks(other_range.start, other_range.stop, from_anchor=False)
        differences = [
            ((len(block.anchor_range) if block.anchor_range else len(block.other_range)), block.category)
            for block in other_blocks
        ]
        proj_feat = ProjectedFeature(
            feature = anchor_feature.feature,
            protein = self.other,
            protein_start = other_range.start + 1,
            protein_stop = other_range.stop,
            reference = False,
            anchor = anchor_feature,
            _differences = ','.join(f'{t[0]}{t[1]}' for t in differences)
        )
        return proj_feat

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
        
        boundaries = [] # Store boundaries of aligned codon blocks
        overhangs = [] # Store placeholders for each codon block, referred to as 'placeholders' in the manuscript.
        categories = [] # Store categories of codon block alignments
        # Define a dictionary to map frame shifts to codon alignment categories
        frame_to_category = {
            0: CodonAlignCat.MATCH,
            1: CodonAlignCat.FRAME_AHEAD,
            2: CodonAlignCat.FRAME_BEHIND
        }
        while tx_blocks:
            tx_category, anchor_tx_range, other_tx_range = tx_blocks.popleft()

            # Check if the block overlaps (derived from transcript blocks) an ORF boundary and split it up if necessary
            if anchor_tx_range.start < 0 < anchor_tx_range.stop:
                # Split the anchor (of reference isoform) range and corresponding other range at the ORF boundary
                anchor_tx_range, next_anchor_range, other_tx_range, next_other_range = split_paired_ranges(anchor_tx_range, other_tx_range, 0)
                tx_blocks.appendleft((tx_category, next_anchor_range, next_other_range))
            if other_tx_range.start < 0 < other_tx_range.stop:
                # Split the other (of alternate isoform) range and corresponding anchor range at the ORF boundary
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
            if (outside_anchor_orf and outside_other_orf
            or outside_anchor_orf and tx_category is SeqAlignCat.DELETION
            or outside_other_orf and tx_category is SeqAlignCat.INSERTION):
                continue
            
            # Convert block range to protein coordinates
            if outside_anchor_orf < 0:
                anchor_pr_start, anchor_start_overhang, anchor_pr_stop, anchor_stop_overhang = 0, 0, 0, 0
            elif outside_anchor_orf > 0:
                anchor_pr_start, anchor_start_overhang, anchor_pr_stop, anchor_stop_overhang = (anchor_pr_len, 0, anchor_pr_len, 0)
            else:
                anchor_pr_start, anchor_start_overhang = divmod(anchor_tx_range.start, 3)
                anchor_pr_stop, anchor_stop_overhang = divmod(anchor_tx_range.stop, 3)
            if outside_other_orf < 0:
                other_pr_start, other_start_overhang, other_pr_stop, other_stop_overhang = 0, 0, 0, 0
            elif outside_other_orf > 0:
                other_pr_start, other_start_overhang, other_pr_stop, other_stop_overhang = (other_pr_len, 0, other_pr_len, 0)
            else:
                other_pr_start, other_start_overhang = divmod(other_tx_range.start, 3)
                other_pr_stop, other_stop_overhang = divmod(other_tx_range.stop, 3)
            # Similar calculations for other_pr_start, other_start_placeholder, other_pr_stop, other_stop_placeholder
          
            # infer codon block category
            if tx_category is SeqAlignCat.MATCH:
                if outside_anchor_orf:
                    cd_category = CodonAlignCat.TRANSLATED
                elif outside_other_orf:
                    cd_category = CodonAlignCat.UNTRANSLATED
                else:
                    frameshift = (other_start_overhang - anchor_start_overhang) % 3
                    cd_category = frame_to_category[frameshift]
            elif tx_category is SeqAlignCat.DELETION:
                cd_category = CodonAlignCat.DELETION
            elif tx_category is SeqAlignCat.INSERTION:
                cd_category = CodonAlignCat.INSERTION
            else:
                raise RuntimeError
              
            # Store the boundaries, placeholders, and categories for each codon block
            boundaries.append((anchor_pr_stop, other_pr_stop))
            overhangs.append((anchor_stop_overhang, other_stop_overhang))
            categories.append(cd_category)
        
        # second pass to adjust edges
        assert overhangs[-1] == (0, 0)
        anchor_boundary_shifts, other_boundary_shifts = dict(), dict()

        i = 0
        while i < len(boundaries) - 1:
            curr_category, next_category = categories[i:i+2]
            overhang = overhangs[i]

            try:
                anchor_boundary = anchor_boundary_shifts[boundaries[i][0]]
            except KeyError:
                anchor_boundary = boundaries[i][0]
                shift_anchor_boundary = (
                    overhang[0] == 2 and (
                        next_category in {CodonAlignCat.MATCH, CodonAlignCat.FRAME_AHEAD}
                        or curr_category is CodonAlignCat.FRAME_AHEAD
                        or curr_category in {CodonAlignCat.DELETION, CodonAlignCat.UNTRANSLATED} and next_category is CodonAlignCat.INSERTION
                    ) or overhang[0] == 1 and curr_category is CodonAlignCat.DELETION and next_category is CodonAlignCat.FRAME_AHEAD
                )
                if shift_anchor_boundary:
                    anchor_boundary_shifts[anchor_boundary] = anchor_boundary + 1
                    anchor_boundary += 1
            try:
                other_boundary = other_boundary_shifts[boundaries[i][1]]
            except KeyError:
                other_boundary = boundaries[i][1]
                shift_other_boundary = (
                    overhang[1] == 2 and (
                        next_category in {CodonAlignCat.MATCH, CodonAlignCat.FRAME_BEHIND}
                        or curr_category is CodonAlignCat.FRAME_BEHIND
                        or curr_category in {CodonAlignCat.INSERTION, CodonAlignCat.TRANSLATED} and next_category is CodonAlignCat.DELETION
                    ) or overhang[1] == 1 and curr_category is CodonAlignCat.INSERTION and next_category is CodonAlignCat.FRAME_BEHIND
                )
                if shift_other_boundary:
                    other_boundary_shifts[other_boundary] = other_boundary + 1
                    other_boundary += 1
            boundaries[i] = anchor_boundary, other_boundary

            # insert a single-codon block if necessary
            category_to_insert = None
            if (curr_category is CodonAlignCat.MATCH and overhang == (2, 2) or next_category is CodonAlignCat.MATCH and overhang == (1, 1)):
                category_to_insert = CodonAlignCat.EDGE
            if (curr_category is CodonAlignCat.FRAME_AHEAD and next_category is CodonAlignCat.DELETION and overhang == (1, 2)
            or curr_category is CodonAlignCat.INSERTION and next_category is CodonAlignCat.FRAME_AHEAD and overhang == (1, 2)
            or curr_category is CodonAlignCat.FRAME_BEHIND and next_category is CodonAlignCat.INSERTION and overhang == (2, 1)
            or curr_category is CodonAlignCat.DELETION and next_category is CodonAlignCat.FRAME_BEHIND and overhang == (2, 1)):
                category_to_insert = CodonAlignCat.COMPLEX

            if category_to_insert:
                boundaries.insert(i+1, (anchor_boundary + 1, other_boundary + 1))
                overhangs.insert(i+1, overhang)
                categories.insert(i+1, category_to_insert)
                anchor_boundary_shifts[anchor_boundary] = anchor_boundary + 1
                other_boundary_shifts[other_boundary] = other_boundary + 1
                i += 1
            i += 1
        #endregion

        # merge consecutive codon blocks w/ same category
        assert len(boundaries) == len(categories)
        cd_blocks: list['CodonAlignmentBlock'] = []
        prev_boundary = (0, 0)
        for category, group in groupby(range(len(categories)), key=categories.__getitem__):
            anchor_start, other_start = prev_boundary
            idx = last(group)
            anchor_stop, other_stop = boundaries[idx]
            prev_boundary = anchor_stop, other_stop
            cblock = CodonAlignmentBlock(range(anchor_start, anchor_stop), range(other_start, other_stop), category=category)
            cd_blocks.append(cblock)
        
        anchor_blocks = IntervalTree.from_tuples(
            (block.anchor_range.start, block.anchor_range.stop, block)
            for block in cd_blocks if block.anchor_range
        )
        other_blocks = IntervalTree.from_tuples(
            (block.other_range.start, block.other_range.stop, block)
            for block in cd_blocks if block.other_range
        )

        # map each cblock to one or more tblocks
        tblock_to_cblocks: dict['TranscriptAlignmentBlock', tuple['CodonAlignmentBlock', ...]] = dict()
        cblock_to_tblock: dict['CodonAlignmentBlock', 'TranscriptAlignmentBlock'] = dict()
        
        def link_cblock_and_tblock(cblock, tblock):
            cblock_to_tblock[cblock] = tblock
            tblock_to_cblocks.setdefault(tblock, []).append(cblock)

        nt_shared, nt_exclusive = partition(lambda t: t[1].category in CD_DEL_INS, enumerate(cd_blocks))
        non_frameshift, frameshift = partition(lambda t: t[1].category in FRAMESHIFT, nt_shared)
        for i, cblock in nt_exclusive:
            if cblock.anchor_range:
                tx_start, tx_stop = map(anchor.get_transcript_coord_from_protein_coord, (cblock.anchor_range.start, cblock.anchor_range.stop))
                mapper = tx_aln.anchor_blocks
            else:
                tx_start, tx_stop = map(other.get_transcript_coord_from_protein_coord, (cblock.other_range.start, cblock.other_range.stop))
                mapper = tx_aln.other_blocks
            intervals = mapper.overlap(tx_start + 1, tx_stop + 1)  # use coord of the nucleotide in the middle of the codon to avoid edge cases
            tblock = one(interval.data for interval in intervals if interval.data.category in SEQ_DEL_INS)
            link_cblock_and_tblock(cblock, tblock)
        
        # map frameshifts to closest preceding del/ins tblock with length indivisible by 3
        match_tblock_to_nonsymmetric_tblock = dict()
        latest_nonsymmetric_tblock = None
        for tblock in tx_aln.blocks:
            if tblock.category is SeqAlignCat.MATCH:
                if latest_nonsymmetric_tblock:
                    match_tblock_to_nonsymmetric_tblock[tblock] = latest_nonsymmetric_tblock
            else:
                if (len(tblock.other_range) - len(tblock.anchor_range)) % 3:
                    latest_nonsymmetric_tblock = tblock
        for i, cblock in frameshift:
            tx_start, tx_stop = map(anchor.get_transcript_coord_from_protein_coord, (cblock.anchor_range.start, cblock.anchor_range.stop))
            intervals = tx_aln.anchor_blocks.overlap(tx_start + 1, tx_stop + 1)
            match_tblock = one(interval.data for interval in intervals if interval.data.category is SeqAlignCat.MATCH)
            try:
                tblock = match_tblock_to_nonsymmetric_tblock[match_tblock]
            except KeyError:
                pass
            else:
                link_cblock_and_tblock(cblock, tblock)
        
        for i, cblock in non_frameshift:
            tblock = None
            if cblock.category in {CodonAlignCat.UNTRANSLATED, CodonAlignCat.TRANSLATED}:
                if not cblock.other_range:
                    nterminal = cblock.other_range.start == 0
                else:
                    nterminal = cblock.anchor_range.start == 0
                if nterminal:
                    mapper = anchor_blocks if cblock.category is CodonAlignCat.UNTRANSLATED else other_blocks
                    start_codon_cblock = only(
                        interval.data for interval in mapper.at(0)
                        if interval.data.category in CD_DEL_INS
                    )
                    tblock = cblock_to_tblock.get(start_codon_cblock)
                else:
                    if cblock.category is CodonAlignCat.UNTRANSLATED:
                        mapper = other_blocks
                        stop_codon_pos = other.length - 1
                    else:
                        mapper = anchor_blocks
                        stop_codon_pos = anchor.length - 1
                    stop_codon_cblock = one(interval.data for interval in mapper.at(stop_codon_pos))
                    tblock = cblock_to_tblock.get(stop_codon_cblock)
            elif cblock.category in SPLIT_CODON:
                del_ins_cblock = cd_blocks[i-1] if cd_blocks[i-1].category in CD_DEL_INS else cd_blocks[i+1]
                tblock = cblock_to_tblock.get(del_ins_cblock)
            if tblock:
                link_cblock_and_tblock(cblock, tblock)

        tblock_to_cblocks = dict((k, tuple(v)) for (k, v) in sorted(tblock_to_cblocks.items()))
        cblock_to_tblock = dict(sorted(cblock_to_tblock.items()))

        return cls(anchor, other, anchor_blocks, other_blocks, tblock_to_cblocks, cblock_to_tblock)


@define(eq=False)
class ProteinAlignment:
    anchor: 'Protein'
    other: 'Protein'
    anchor_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False, validator=check_block_ranges)
    other_blocks: 'IntervalTree' = field(factory=IntervalTree, repr=False, validator=check_block_ranges)
    cblock_to_pblock: dict['CodonAlignmentBlock', 'ProteinAlignmentBlock'] = field(factory=dict, repr=False)
    pblock_to_cblocks: dict['ProteinAlignmentBlock', tuple['CodonAlignmentBlock', ...]] = field(factory=dict, repr=False)

    def __attrs_post_init__(self):
        for pblock in self.blocks:
            cblocks = self.pblock_to_cblocks.get(pblock, ())
            if any(pblock is not self.cblock_to_pblock.get(cblock) for cblock in cblocks):
                raise ValueError(f'Broken pblock-to-cblock mapping for {pblock}')

    @property
    def blocks(self) -> tuple['ProteinAlignmentBlock', ...]:
        return tuple(sorted({i.data for i in chain(self.anchor_blocks, self.other_blocks)}))

    @classmethod
    @lru_cache(maxsize=CACHE_SIZE)
    def from_proteins(cls, anchor: 'Protein', other: 'Protein'):
        cd_aln = CodonAlignment.from_proteins(anchor, other)

        pblock_to_cblocks: dict['ProteinAlignmentBlock', tuple['CodonAlignmentBlock', ...]] = dict()
        cblock_to_pblock: dict['CodonAlignmentBlock', 'ProteinAlignmentBlock'] = dict()
        
        pblock_bounds: list[int] = [
            last(cblock_indices) + 1
            for _, cblock_indices in groupby(
                range(len(cd_aln.blocks)),
                key = lambda i: cd_aln.blocks[i].category is CodonAlignCat.MATCH
            )
        ]
        pblock_categories = []
        pblock_split5 = [False for pblock in pblock_bounds]
        pblock_split3 = [False for pblock in pblock_bounds]
        pblock_ragged5 = [False for pblock in pblock_bounds]
        pblock_ragged3 = [False for pblock in pblock_bounds]
        for p, (c0, c1) in enumerate(windowed([0] + pblock_bounds, 2)):
            cblocks = cd_aln.blocks[c0:c1]
            anchor_start, anchor_stop = cblocks[0].anchor_range.start, cblocks[-1].anchor_range.stop
            other_start, other_stop = cblocks[0].other_range.start, cblocks[-1].other_range.stop
            reduced_cblock_categories = {cblock.category for cblock in cblocks} - SPLIT_CODON
            anchor_sequence = anchor.sequence[anchor_start:anchor_stop]
            other_sequence = other.sequence[other_start:other_stop]
            if anchor_sequence == other_sequence:
                category = SeqAlignCat.MATCH
            elif reduced_cblock_categories <= ANCHOR_EXCLUSIVE:
                category = SeqAlignCat.DELETION
                pblock_split5[p] = cblocks[0].category in SPLIT_CODON
                pblock_split3[p] = cblocks[-1].category in SPLIT_CODON
            elif reduced_cblock_categories <= OTHER_EXCLUSIVE:
                category = SeqAlignCat.INSERTION
                pblock_split5[p] = cblocks[0].category in SPLIT_CODON
                pblock_split3[p] = cblocks[-1].category in SPLIT_CODON
            else:
                category = SeqAlignCat.SUBSTITUTION
            pblock_categories.append(category)
            pblock_ragged5[p] = pblock_split5[p] and anchor_sequence[0] != other_sequence[0]
            pblock_ragged3[p] = pblock_split3[p] and anchor_sequence[-1] != other_sequence[-1]
        
        assert len(pblock_bounds) == len(pblock_categories)

        # second pass to move synonymous split codons into match pblocks
        for p, (split5, split3, ragged5, ragged3) in enumerate(zip(pblock_split5, pblock_split3, pblock_ragged5, pblock_ragged3)):
            if split5 and not ragged5:
                pblock_bounds[p-1] += 1
            if split3 and not ragged3:
                pblock_bounds[p] -= 1

        for p, (c0, c1) in enumerate(windowed([0] + pblock_bounds, 2)):
            cblocks = cd_aln.blocks[c0:c1]
            anchor_range = range(cblocks[0].anchor_range.start, cblocks[-1].anchor_range.stop)
            other_range = range(cblocks[0].other_range.start, cblocks[-1].other_range.stop)
            pblock = ProteinAlignmentBlock(
                anchor_range,
                other_range,
                category = pblock_categories[p],
                ragged5 = pblock_ragged5[p],
                ragged3 = pblock_ragged3[p]
            )
            pblock_to_cblocks[pblock] = cblocks
            for cblock in cblocks:
                cblock_to_pblock[cblock] = pblock
        # TODO: second pass to merge match pblocks

        anchor_blocks = IntervalTree.from_tuples((block.anchor_range.start, block.anchor_range.stop, block) for block in pblock_to_cblocks if block.anchor_range)
        other_blocks = IntervalTree.from_tuples((block.other_range.start, block.other_range.stop, block) for block in pblock_to_cblocks if block.other_range)

        return cls(anchor, other, anchor_blocks, other_blocks, cblock_to_pblock, pblock_to_cblocks)
