from abc import ABC
from collections import deque
from collections.abc import Sequence
from itertools import chain, groupby
from operator import attrgetter
from typing import Iterable, List, Optional, Union, MutableSequence

from constants import AminoAcid, ProteinLevelAlignmentCategory
from constants import TranscriptLevelAlignmentCategory as TranscriptAlignCat
from constants import ProteinLevelAlignmentCategory as ProteinAlignCat
from models import ORF, Exon, Nucleotide, Protein, Residue, Strand, Transcript


def get_first_nt_adjusted_coord(res: 'Residue', strand: 'Strand' = Strand.PLUS) -> int:
            if len(res.exons) == 1 or res.codon[0].exon is res.codon[1].exon:
                return res.codon[0].coordinate
            # if codon's primary exon is downstream, "abacus" the coord of the 1st nucleotide
            elif strand is Strand.PLUS:
                return res.codon[1].coordinate - 1
            elif strand is Strand.MINUS:
                return res.codon[1].coordinate + 1


class GapResidue(Residue):
    def __init__(self, protein: 'Protein', position: int, upstream_exon: 'Exon', downstream_exon: 'Exon'):
        # if upstream_exon.transcript is not downstream_exon.transcript:
        #     raise ValueError(f'{upstream_exon} and {downstream_exon} belong to different transcripts')
        # if upstream_exon.transcript is not protein.orf.transcript:
        #     raise ValueError(f'{upstream_exon} and {protein} belong to different transcripts')
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
    def __init__(self, parent, position, start, end, category):
        self.parent = parent
        self.position = position
        self.start = start  # 0-based, inclusive
        self.end = end  # 0-based, exclusive
        self.length = end - start
        self.category = category
        self.annotation = None

    def __repr__(self):
        return f'{repr(self.parent)}:{self.position}-{self.category}'

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
        ttype_str = ''.join(str(res.category) for res in self)
        return anchor_str + '\n' + ttype_str + '\n' + other_str
    

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
        self.protein_blocks = get_protein_blocks(self, self.transcript_blocks)

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
    
    def annotate(self):
        pass
    

def rough_alignment(anchor: 'Protein', other: 'Protein', strand: 'Strand') -> List['ResidueAlignment']:
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
                elif abs(coord_diff) == 1:
                    event_type = TranscriptAlignCat.FRAMESHIFT
                elif coord_diff < 0:
                    event_type = TranscriptAlignCat.DELETION
                elif coord_diff > 0:
                    event_type = TranscriptAlignCat.INSERTION
            else:
                # for split codons, need to compare the coords of all their nucleotides
                anchor_coords = tuple(nt.coordinate for nt in anchor_current.codon)
                other_coords = tuple(nt.coordinate for nt in other_current.codon)
                overlap = len(set(anchor_coords) & set(other_coords))
                if overlap < 2:
                    coord_diff = anchor_coords[1] - other_coords[1]
                    if strand is Strand.MINUS:
                        coord_diff = -coord_diff
                    if coord_diff < 0:
                        event_type = TranscriptAlignCat.DELETION
                    elif coord_diff > 0:
                        event_type = TranscriptAlignCat.INSERTION
                elif overlap == 2:
                    if anchor_coords[1] != other_coords[1]:
                        event_type = TranscriptAlignCat.FRAMESHIFT
                    elif anchor_current.amino_acid is not other_current.amino_acid:
                        event_type = TranscriptAlignCat.EDGE_MISMATCH
                elif overlap == 3 and anchor_current.amino_acid is other_current.amino_acid:
                    event_type = TranscriptAlignCat.MATCH
        elif not other_stack:
            event_type = TranscriptAlignCat.DELETION
        elif not anchor_stack:
            event_type = TranscriptAlignCat.INSERTION
        
        def make_gap_res_from_prev_res(prev_res: 'Residue', protein: 'Protein'):
            gap_pos = prev_res.position if prev_res else 0
            upstream_exon = prev_res.codon[0].exon if prev_res else None
            return GapResidue(protein, gap_pos, upstream_exon, None)
        def get_next_res_from_stack(stack: MutableSequence['Residue'], gap_res: 'GapResidue'):
            next_res = stack.popleft()
            if gap_res:
                gap_res.downstream_exon = next_res.codon[-1].exon
            return next_res

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
    delete_insert = {TranscriptAlignCat.DELETION, TranscriptAlignCat.INSERTION}
    i = 1
    while i+1 < len(chain):
        curr = chain[i]
        curr_cat = curr.category
        if curr_cat in delete_insert:
            prev = chain[i-1]
            next = chain[i+1]
            prev_cat = prev.category
            next_cat = next.category
            if len({prev_cat, curr_cat, next_cat}) == 3:
                merge_res = None
                if prev_cat in delete_insert and next_cat is TranscriptAlignCat.FRAMESHIFT:
                    merge_res = prev
                elif next_cat in delete_insert and prev_cat is TranscriptAlignCat.FRAMESHIFT:
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


def get_transcript_blocks(aln: 'TranscriptBasedAlignment') -> List['AlignmentBlock']:
    tblocks = []
    start = 0
    for i, (category, res_alns) in enumerate(groupby(aln, key=attrgetter('category'))):
        tblock_length = len(list(res_alns))
        tblocks.append(AlignmentBlock(aln, i, start, start+tblock_length, category))
        start += tblock_length
    return tblocks


def get_protein_blocks(parent, tblocks):
    # TODO: account for amino acid sequence
    pblocks = []
    for i, (is_match, tblock_group) in enumerate(groupby(tblocks, key=lambda tblock: tblock.category is TranscriptAlignCat.MATCH), start=1):
        tblock_group = list(tblock_group)
        start = min(tblock.start for tblock in tblock_group)
        end = max(tblock.end for tblock in tblock_group)
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
        pblocks.append(AlignmentBlock(parent, i, start, end, pblock_category))
    return pblocks
    