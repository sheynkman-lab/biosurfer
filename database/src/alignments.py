from abc import ABC
from collections import deque
from typing import List, Optional, Union

from constants import AminoAcid
from constants import TranscriptLevelAlignmentCategory as TranscriptAlignCat
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


class TranscriptBasedAlignment(Alignment):
    def __init__(self, anchor: 'Protein', other: 'Protein'):
        if anchor.orf.gene is not other.orf.gene:
            raise ValueError(f'{anchor} and {other} belong to different genes')
        if anchor.orf.transcript.strand is not other.orf.transcript.strand:
            raise ValueError(f'{anchor.orf.transcript} and {other.orf.transcript} are on different strands')
        else:
            strand = anchor.orf.transcript.strand
        super().__init__(anchor, other)
        
        anchor_stack = deque(anchor.residues)
        other_stack = deque(other.residues)
        anchor_res, other_res = None, None
        anchor_gap, other_gap = None, None
        self.chain: List['ResidueAlignment'] = []
        while anchor_stack or other_stack:
            # use residues at top of stack to determine event type
            event_type = TranscriptAlignCat.OTHER
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

            # TODO: set upstream and downstream exons for GapResidue
            if event_type is TranscriptAlignCat.DELETION:
                anchor_gap = None
                if not other_gap:
                    gap_pos = other_res.position if other_res else 0
                    other_gap = GapResidue(other, gap_pos, None, None)
                anchor_res = anchor_stack.popleft()
                res_align = ResidueAlignment(anchor_res, other_gap, TranscriptAlignCat.DELETION)
            elif event_type is TranscriptAlignCat.INSERTION:
                other_gap = None
                if not anchor_gap:
                    gap_pos = anchor_res.position if anchor_res else 0
                    anchor_gap = GapResidue(anchor, gap_pos, None, None)
                other_res = other_stack.popleft()
                res_align = ResidueAlignment(anchor_gap, other_res, TranscriptAlignCat.INSERTION)
            else:
                anchor_gap, other_gap = None, None
                anchor_res = anchor_stack.popleft()
                other_res = other_stack.popleft()
                res_align = ResidueAlignment(anchor_res, other_res, event_type)
            
            assert res_align.anchor.protein is anchor
            assert res_align.other.protein is other
            self.chain.append(res_align)
    
    @property
    def full(self):
        anchor_str = ''.join(str(res.anchor.amino_acid) for res in self.chain)
        other_str = ''.join(str(res.other.amino_acid) for res in self.chain)
        ttype_str = ''.join(str(res.category) for res in self.chain)
        return anchor_str + '\n' + ttype_str + '\n' + other_str
    


    
