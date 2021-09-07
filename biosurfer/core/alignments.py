from abc import ABC, abstractmethod
from collections import deque
from collections.abc import Sequence
from itertools import chain, groupby
from operator import attrgetter
from typing import Iterable, List, Optional, Union, MutableSequence

from .constants import AminoAcid, ProteinRegion, Strand
from .constants import TranscriptLevelAlignmentCategory as TranscriptAlignCat
from .constants import ProteinLevelAlignmentCategory as ProteinAlignCat
from .models import Transcript, Exon, Nucleotide, Protein, Residue, Transcript


# def get_first_nt_adjusted_coord(res: 'Residue', strand: 'Strand' = Strand.PLUS) -> int:
#             if len(res.exons) == 1 or res.codon[0].exon is res.codon[1].exon:
#                 return res.codon[0].coordinate
#             # if codon's primary exon is downstream, "abacus" the coord of the 1st nucleotide
#             elif strand is Strand.PLUS:
#                 return res.codon[1].coordinate - 1
#             elif strand is Strand.MINUS:
#                 return res.codon[1].coordinate + 1


class GapResidue(Residue):
    def __init__(self, protein: 'Protein', position: int, upstream_exon: Optional['Exon'], downstream_exon: Optional['Exon']):
        super().__init__(protein, AminoAcid.GAP, position)
        self.upstream_exon = upstream_exon
        self.downstream_exon = downstream_exon
    
    @Residue.codon_str.getter
    def codon_str(self):
        return '-'

    @Residue.exons.getter
    def exons(self):
        return [self.upstream_exon, self.downstream_exon]
    
    @Residue.primary_exon.getter
    def primary_exon(self):
        return None


# TODO: need to clean up Alignment class hierarchy
class Alignment(ABC):
    def __init__(self, anchor, other):
        self.anchor = anchor
        self.other = other
    
    def __repr__(self):
        return f'{self.anchor}|{self.other}'


class ResidueAlignment(Alignment):
    def __init__(self, anchor: 'Residue', other: 'Residue', category: 'TranscriptAlignCat'):
        super().__init__(anchor, other)
        self.category = category


class AlignmentBlock(Sequence):
    def __init__(self, parent, position, start, end):
        self.parent = parent
        self.position = position
        self.start = start  # 0-based, inclusive
        self.end = end  # 0-based, exclusive
        self.length = end - start

    @abstractmethod
    def __repr__(self):
        pass

    def __len__(self):
        return self.length

    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.step:
                raise NotImplementedError('AlignmentBlock does not support slices with steps')
            return [self[i] for i in range(*index.indices(self.length))]
        elif isinstance(index, int):
            return self.parent[self._parent_based_index(index)]
        else:
            raise TypeError(f'{index} is not an int or slice')
    
    def _parent_based_index(self, index):
        if index < -self.length or index >= self.length:
            raise IndexError(f'{index} out of range for {self}')
        return (index % self.length) + self.start
    
    @property
    def full(self):
        anchor_str = ''.join(str(res.anchor.amino_acid) for res in self)
        other_str = ''.join(str(res.other.amino_acid) for res in self)
        cat_str = ''.join(str(res.category) for res in self)
        return anchor_str + '\n' + cat_str + '\n' + other_str
    

class TranscriptAlignmentBlock(AlignmentBlock):
    def __init__(self, parent, position, start, end, category: 'TranscriptAlignCat'):
        super().__init__(parent, position, start, end)
        self.category = category
        # These attributes are useful in the annotation code
        self._prev_match_or_frame_tblock = None
        self._next_match_or_frame_tblock = None
    
    def __repr__(self):
        return f'{self.parent}:tblock{self.position}-{self.category}'


class ProteinAlignmentBlock(AlignmentBlock):
    def __init__(self, parent, position, tblocks: Sequence['TranscriptAlignmentBlock'], category: 'ProteinAlignCat'):
        start = min(tblock.start for tblock in tblocks)
        end = max(tblock.end for tblock in tblocks)
        super().__init__(parent, position, start, end)
        self.category = category
        self.transcript_blocks = list(tblocks)
        self._annotations = []
        self.region = ProteinRegion.INTERNAL

    def __repr__(self):
        return f'{self.parent}:pblock{self.position}-{self.category}'
    
    @property
    def annotation(self):
        return ', \n'.join(self._annotations) if self._annotations else None


class TranscriptBasedAlignment(Alignment, Sequence):
    def __init__(self, anchor: 'Protein', other: 'Protein'):
        if anchor.orf.gene is not other.orf.gene:
            raise ValueError(f'{anchor} and {other} belong to different genes')
        if anchor.orf.transcript.strand is not other.orf.transcript.strand:
            raise ValueError(f'{anchor.orf.transcript} and {other.orf.transcript} are on different strands')
        else:
            strand = anchor.orf.transcript.strand
        super().__init__(anchor, other)
        self._chain = rough_alignment(anchor, other, strand)
        refine_alignment(self._chain)
        self.transcript_blocks = get_transcript_blocks(self)
        self.protein_blocks = get_protein_blocks(self)

    def __getitem__(self, index):
        return self._chain[index]
    
    def __len__(self):
        return len(self._chain)

    @property
    def full(self):
        anchor_str = ''.join(str(res.anchor.amino_acid) for res in self)
        other_str = ''.join(str(res.other.amino_acid) for res in self)
        ttype_str = ''.join(str(res.category) for res in self)
        return anchor_str + '\n' + ttype_str + '\n' + other_str
    
    # TODO: consider using Annotation classes in the future?
    # TODO: detect NAGNAG splicing
    def annotate(self) -> None:
        FRAMESHIFT = {TranscriptAlignCat.FRAME_AHEAD, TranscriptAlignCat.FRAME_BEHIND}
        # DELETE_INSERT = {TranscriptAlignCat.DELETION, TranscriptAlignCat.INSERTION}

        strand = self.anchor.orf.transcript.strand
        anchor_transcript = self.anchor.orf.transcript
        other_transcript = self.other.orf.transcript
        nterminal_pblock = None
        cterminal_pblock = None
        upstream_cterm_tblock = None
        upstream_cterm_res_aln = None
        
        # classify internal splicing events
        for pblock in self.protein_blocks:
            if not nterminal_pblock:
                for res_aln in pblock:
                    if res_aln.anchor.position == 1 or res_aln.other.position == 1:
                        nterminal_pblock = pblock
                        pblock.region = ProteinRegion.NTERMINUS
                        break
            if not cterminal_pblock:
                for res_aln in pblock:
                    if res_aln.anchor.amino_acid is AminoAcid.STOP or res_aln.other.amino_acid is AminoAcid.STOP:
                        cterminal_pblock = pblock
                        pblock.region = ProteinRegion.CTERMINUS
                        for tblock in pblock.transcript_blocks:
                            if res_aln in tblock:
                                upstream_cterm_tblock = tblock
                                break
                        upstream_cterm_res_aln = res_aln
                        break
            
            if pblock.category is ProteinAlignCat.MATCH:
                continue
            
            for tblock in pblock.transcript_blocks:
                if tblock.category is TranscriptAlignCat.DELETION:
                    first_exon = tblock[0].anchor.codon[2].exon
                    last_exon = tblock[-1].anchor.codon[0].exon

                    if tblock._prev_match_or_frame_tblock and tblock._next_match_or_frame_tblock:
                        prev_anchor_exon = tblock._prev_match_or_frame_tblock[-1].anchor.codon[0].exon
                        next_anchor_exon = tblock._next_match_or_frame_tblock[0].anchor.codon[2].exon
                        if prev_anchor_exon is next_anchor_exon:
                            pblock._annotations.append(f'portion of {prev_anchor_exon} intronized')
                        else:
                            e_first = first_exon.position
                            e_last = last_exon.position
                            if prev_anchor_exon is first_exon:
                                pblock._annotations.append(f'{first_exon} shortened by alternative splice donor')
                                e_first += 1
                            if next_anchor_exon is last_exon:
                                pblock._annotations.append(f'{last_exon} shortened by alternative splice acceptor')
                                e_last -= 1
                            first_skipped_exon = anchor_transcript.exons[e_first-1]
                            last_skipped_exon = anchor_transcript.exons[e_last-1]
                            if first_skipped_exon is last_skipped_exon:
                                pblock._annotations.append(f'{first_skipped_exon} skipped')
                            elif e_first < e_last:
                                pblock._annotations.append(f'exons {first_skipped_exon} to {last_skipped_exon} skipped')
                
                elif tblock.category is TranscriptAlignCat.INSERTION:
                    first_exon = tblock[0].other.codon[2].exon
                    last_exon = tblock[-1].other.codon[0].exon

                    if tblock._prev_match_or_frame_tblock and tblock._next_match_or_frame_tblock:
                        prev_anchor_exon = tblock._prev_match_or_frame_tblock[-1].anchor.codon[0].exon
                        next_anchor_exon = tblock._next_match_or_frame_tblock[0].anchor.codon[2].exon
                        prev_other_exon = tblock._prev_match_or_frame_tblock[-1].other.codon[0].exon
                        next_other_exon = tblock._next_match_or_frame_tblock[0].other.codon[2].exon
                        if prev_other_exon is next_other_exon:
                            pblock._annotations.append(f'retained intron between {prev_anchor_exon} and {next_anchor_exon}')
                        else:
                            number_of_included_exons = last_exon.position - first_exon.position + 1
                            if prev_other_exon is first_exon:
                                pblock._annotations.append(f'{prev_anchor_exon} lengthened by alternative splice donor')
                                number_of_included_exons -= 1
                            if next_other_exon is last_exon:
                                pblock._annotations.append(f'{next_anchor_exon} lengthened by alternative splice acceptor')
                                number_of_included_exons -= 1
                            if number_of_included_exons == 1:
                                pblock._annotations.append(f'exon included between {prev_anchor_exon} and {next_anchor_exon}')
                            elif number_of_included_exons > 1:
                                pblock._annotations.append(f'{number_of_included_exons} exons included between {prev_anchor_exon} and {next_anchor_exon}')
                
                elif tblock.category is TranscriptAlignCat.EDGE_MISMATCH:
                    pblock._annotations.append(f'{tblock[0].anchor} replaced with {tblock[0].other} due to use of alternate junction')
                
                elif tblock.category in FRAMESHIFT:
                    first_exon = tblock[0].anchor.codon[2].exon
                    last_exon = tblock[-1].anchor.codon[0].exon
                    if first_exon is last_exon:
                        pblock._annotations.append(f'{first_exon} translated in different frame')
                    else:
                        pblock._annotations.append(f'{first_exon} to {last_exon} translated in different frame')
        
        # classify N-terminal changes (if any)
        if nterminal_pblock.category is not ProteinAlignCat.MATCH:
            upstream_start_codon = tuple(nt.coordinate for nt in self.anchor.residues[0].codon)
            downstream_start_codon = tuple(nt.coordinate for nt in self.other.residues[0].codon)
            upstream_start_transcript, downstream_start_transcript = anchor_transcript, other_transcript
            if strand is Strand.PLUS:
                alternative_tss = anchor_transcript.start != other_transcript.start
                if upstream_start_codon[1] > downstream_start_codon[1]:
                    upstream_start_codon, downstream_start_codon = downstream_start_codon, upstream_start_codon 
                    upstream_start_transcript, downstream_start_transcript = downstream_start_transcript, upstream_start_transcript
            elif strand is Strand.MINUS:
                alternative_tss = anchor_transcript.stop != other_transcript.stop
                if upstream_start_codon[1] < downstream_start_codon[1]:
                    upstream_start_codon, downstream_start_codon = downstream_start_codon, upstream_start_codon 
                    upstream_start_transcript, downstream_start_transcript = downstream_start_transcript, upstream_start_transcript
            
            upstream_start_codon_shared_nts = sum(downstream_start_transcript.contains_coordinate(coord) for coord in upstream_start_codon)
            downstream_start_codon_shared_nts = sum(upstream_start_transcript.contains_coordinate(coord) for coord in downstream_start_codon)
            if upstream_start_codon_shared_nts == 3:
                if downstream_start_codon_shared_nts == 3:
                    # mutually shared start codons
                    if other_transcript is downstream_start_transcript:
                        nterminal_pblock._annotations.append('usage of downstream alternative TIS')  # TODO: indicate anchor exon
                    else:
                        nterminal_pblock._annotations.append('usage of upstream alternative TIS')  # TODO: indicate anchor exon
                else:
                    # shared upstream start, exclusive downstream start
                    if other_transcript is downstream_start_transcript:
                        nterminal_pblock._annotations.append('usage of downstream TIS revealed by splicing')  # TODO: indicate surrounding anchor exons
                    else:
                        nterminal_pblock._annotations.append('usage of upstream alternative TIS due to removal of anchor TIS by splicing')  # TODO: indicate anchor exon
            else:
                if downstream_start_codon_shared_nts == 3:
                    # exclusive upstream start, shared downstream start
                    if strand is Strand.PLUS:
                        cause = 'alternative TSS' if upstream_start_codon[0] < downstream_start_transcript.start else 'splicing'
                    elif strand is Strand.MINUS:
                        cause = 'alternative TSS' if upstream_start_codon[0] > downstream_start_transcript.stop else 'splicing'
                    if other_transcript is downstream_start_transcript:
                        nterminal_pblock._annotations.append('usage of downstream alternative TIS due to removal of anchor TIS by ' + cause)  # TODO: indicate anchor exon
                    else:
                        nterminal_pblock._annotations.append('usage of upstream TIS revealed by ' + cause)  # TODO: indicate surrounding anchor exons
                else:
                    # mutually exclusive start codons
                    if alternative_tss:
                        nterminal_pblock._annotations.append('alternative TSS leading to mutually exclusive start codons')
                    else:
                        nterminal_pblock._annotations.append('5\' UTR splicing leading to mutually exclusive start codons')
        
        # classify C-terminal changes (if any)
        if cterminal_pblock.category is not ProteinAlignCat.MATCH:
            if upstream_cterm_res_aln.category in FRAMESHIFT:
                if upstream_cterm_res_aln.anchor.amino_acid is AminoAcid.STOP:
                    cterminal_pblock._annotations.append('splicing-induced frameshift leading to usage of downstream stop codon')  # TODO: indicate location of other stop codon
                else:
                    cterminal_pblock._annotations.append('splicing-induced frameshift leading to usage of upstream stop codon')  # TODO: indicate location of other stop codon
            elif upstream_cterm_res_aln.category is TranscriptAlignCat.DELETION:
                if strand is Strand.PLUS:
                    alt_cterm_exons = upstream_cterm_res_aln.anchor.exons[-1].stop < self.other.orf.exons[-1].start
                elif strand is Strand.MINUS:
                    alt_cterm_exons = upstream_cterm_res_aln.anchor.exons[-1].start > self.other.orf.exons[-1].stop
                if alt_cterm_exons:
                    cterminal_pblock._annotations.append('alternative C-terminal exon')
                else:
                    cterminal_pblock._annotations.append('anchor stop codon spliced out leading to usage of downstream stop codon')  # TODO: indicate location of other stop codon
            elif upstream_cterm_res_aln.category is TranscriptAlignCat.INSERTION:
                exon_extension_introduces_stop = upstream_cterm_tblock._prev_match_or_frame_tblock[-1].other.codon[0].exon is upstream_cterm_res_aln.other.codon[2].exon
                if strand is Strand.PLUS:
                    alt_cterm_exons = upstream_cterm_res_aln.other.exons[-1].stop < self.anchor.orf.exons[-1].start
                elif strand is Strand.MINUS:
                    alt_cterm_exons = upstream_cterm_res_aln.other.exons[-1].start > self.anchor.orf.exons[-1].stop
                if exon_extension_introduces_stop:
                    lengthened_exon = upstream_cterm_tblock._prev_match_or_frame_tblock[-1].anchor.codon[0].exon
                    cterminal_pblock._annotations.append(f'upstream stop codon introduced by extension of {lengthened_exon}')
                elif alt_cterm_exons:
                    cterminal_pblock._annotations.append('alternative C-terminal exon')
                else:
                    cterminal_pblock._annotations.append('upstream stop codon introduced by splicing')  # TODO: indicate surrounding anchor exons
            else:
                cterminal_pblock._annotations.append('complex C-terminal event')
            if self.other.orf.nmd:
                cterminal_pblock._annotations.append('NMD candidate')


def rough_alignment(anchor: 'Protein', other: 'Protein', strand: 'Strand') -> List['ResidueAlignment']:
    # helper functions
    def make_gap_res_from_prev_res(prev_res: 'Residue', protein: 'Protein'):
        gap_pos = prev_res.position if prev_res else 0
        upstream_exon = prev_res.codon[0].exon if prev_res else None
        return GapResidue(protein, gap_pos, upstream_exon, None)
    def get_next_res_from_stack(stack: MutableSequence['Residue'], gap_res: 'GapResidue'):
        next_res = stack.popleft()
        if gap_res:
            gap_res.downstream_exon = next_res.codon[-1].exon
        return next_res

    anchor_stack = deque(anchor.residues)
    other_stack = deque(other.residues)
    anchor_res, other_res = None, None
    anchor_gap, other_gap = None, None
    chain = []
    while anchor_stack or other_stack:
        # use residues at top of stack to determine event type
        event_type = TranscriptAlignCat.UNKNOWN
        if anchor_stack and other_stack:
            anchor_current, other_current = anchor_stack[0], other_stack[0]
            if len(anchor_current.exons) == 1 and len(other_current.exons) == 1:
                # for contiguous (nonsplit) codons, can just compare the coords of their middle nucleotides
                coord_diff = anchor_current.codon[1].coordinate - other_current.codon[1].coordinate
                if strand is Strand.MINUS:
                    coord_diff = -coord_diff
                if coord_diff == 0:
                    if anchor_current.amino_acid is other_current.amino_acid:
                        event_type = TranscriptAlignCat.MATCH
                elif coord_diff == -1:
                    event_type = TranscriptAlignCat.FRAME_AHEAD
                elif coord_diff == 1:
                    event_type = TranscriptAlignCat.FRAME_BEHIND
                elif coord_diff < 0:
                    event_type = TranscriptAlignCat.DELETION
                elif coord_diff > 0:
                    event_type = TranscriptAlignCat.INSERTION
            else:
                # for split codons, need to compare the coords of all their nucleotides
                anchor_coords = tuple(nt.coordinate for nt in anchor_current.codon)
                other_coords = tuple(nt.coordinate for nt in other_current.codon)
                overlap = len(set(anchor_coords) & set(other_coords))
                coord_diff = anchor_coords[1] - other_coords[1]
                if strand is Strand.MINUS:
                    coord_diff = -coord_diff
                if overlap < 2:
                    if coord_diff < 0:
                        event_type = TranscriptAlignCat.DELETION
                    elif coord_diff > 0:
                        event_type = TranscriptAlignCat.INSERTION
                elif overlap == 2:
                    if coord_diff < 0:
                        event_type = TranscriptAlignCat.FRAME_AHEAD
                    elif coord_diff > 0:
                        event_type = TranscriptAlignCat.FRAME_BEHIND
                    elif anchor_current.amino_acid is not other_current.amino_acid:
                        event_type = TranscriptAlignCat.EDGE_MISMATCH
                elif overlap == 3 and anchor_current.amino_acid is other_current.amino_acid:
                    event_type = TranscriptAlignCat.MATCH
        elif not other_stack:
            event_type = TranscriptAlignCat.DELETION
        elif not anchor_stack:
            event_type = TranscriptAlignCat.INSERTION

        if event_type is TranscriptAlignCat.DELETION:
            anchor_res = get_next_res_from_stack(anchor_stack, anchor_gap)
            anchor_gap = None
            if not other_gap:
                other_gap = make_gap_res_from_prev_res(other_res, other)
            res_align = ResidueAlignment(anchor_res, other_gap, TranscriptAlignCat.DELETION)
        elif event_type is TranscriptAlignCat.INSERTION:
            other_res = get_next_res_from_stack(other_stack, other_gap)
            other_gap = None
            if not anchor_gap:
                anchor_gap = make_gap_res_from_prev_res(anchor_res, anchor)
            res_align = ResidueAlignment(anchor_gap, other_res, TranscriptAlignCat.INSERTION)
        else:
            anchor_res = get_next_res_from_stack(anchor_stack, anchor_gap)
            other_res = get_next_res_from_stack(other_stack, other_gap)
            anchor_gap, other_gap = None, None
            res_align = ResidueAlignment(anchor_res, other_res, event_type)
        
        assert res_align.anchor.protein is anchor
        assert res_align.other.protein is other
        chain.append(res_align)
    return chain


def refine_alignment(chain: List['ResidueAlignment']) -> None:
    """Performs a second pass on an existing ResidueAlignment chain to identify complex codon alignments."""
    DELETE_INSERT = {TranscriptAlignCat.DELETION, TranscriptAlignCat.INSERTION}
    FRAMESHIFT = {TranscriptAlignCat.FRAME_AHEAD, TranscriptAlignCat.FRAME_BEHIND}
    i = 1
    while i+1 < len(chain):
        curr = chain[i]
        curr_cat = curr.category
        if curr_cat in DELETE_INSERT:
            prev = chain[i-1]
            next = chain[i+1]
            prev_cat = prev.category
            next_cat = next.category
            if len({prev_cat, curr_cat, next_cat}) == 3:
                merge_res = None
                if prev_cat in DELETE_INSERT and next_cat in FRAMESHIFT:
                    merge_res = prev
                elif next_cat in DELETE_INSERT and prev_cat in FRAMESHIFT:
                    merge_res = next
                if merge_res:
                    if curr_cat is TranscriptAlignCat.DELETION:
                        chain[i] = ResidueAlignment(curr.anchor, merge_res.other, TranscriptAlignCat.COMPLEX)
                    elif curr_cat is TranscriptAlignCat.INSERTION:
                        chain[i] = ResidueAlignment(merge_res.anchor, curr.other, TranscriptAlignCat.COMPLEX)
                    if merge_res is prev:
                        chain.pop(i-1)
                        continue  # popping previous item from chain automatically advances i to the next item
                    elif merge_res is next:
                        chain.pop(i+1)
        i += 1


def get_transcript_blocks(aln: Iterable['ResidueAlignment']) -> List['TranscriptAlignmentBlock']:
    tblocks = []
    start = 0
    prev_match_or_frame_tblock = None
    MATCH_OR_FRAME = {TranscriptAlignCat.MATCH, TranscriptAlignCat.FRAME_AHEAD, TranscriptAlignCat.FRAME_BEHIND}
    for i, (category, res_alns) in enumerate(groupby(aln, key=attrgetter('category'))):
        tblock_length = len(list(res_alns))
        tblock = TranscriptAlignmentBlock(aln, i, start, start+tblock_length, category)
        if category in MATCH_OR_FRAME:
            for prev_tblock in reversed(tblocks):
                if prev_tblock is prev_match_or_frame_tblock:
                    break
                prev_tblock._next_match_or_frame_tblock = tblock
            prev_match_or_frame_tblock = tblock
        else:
            tblock._prev_match_or_frame_tblock = prev_match_or_frame_tblock
        tblocks.append(tblock)
        start += tblock_length
    return tblocks


def get_protein_blocks(parent: 'TranscriptBasedAlignment') -> List['ProteinAlignmentBlock']:
    # TODO: account for amino acid sequence
    pblocks = []
    for i, (is_match, tblock_group) in enumerate(groupby(parent.transcript_blocks, key=lambda tblock: tblock.category is TranscriptAlignCat.MATCH)):
        tblock_group = list(tblock_group)
        if is_match:
            pblock_category = ProteinAlignCat.MATCH
        else:
            categories = {tblock.category for tblock in tblock_group if tblock.category is not TranscriptAlignCat.EDGE_MISMATCH}
            pblock_category = ProteinAlignCat.SUBSTITUTION
            if len(categories) == 1:
                single_category = list(categories)[0]
                if single_category is TranscriptAlignCat.DELETION:
                    pblock_category = ProteinAlignCat.DELETION
                elif single_category is TranscriptAlignCat.INSERTION:
                    pblock_category = ProteinAlignCat.INSERTION
        pblocks.append(ProteinAlignmentBlock(parent, i, tblock_group, pblock_category))
    return pblocks
    