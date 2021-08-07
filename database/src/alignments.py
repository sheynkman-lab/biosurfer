from abc import ABC
from collections import deque
from math import inf, isfinite
from typing import Optional

from constants import AminoAcid, ProteinLevelEvent, TranscriptLevelEvent
from models import ORF, Exon, Nucleotide, Protein, Residue, Strand, Transcript


def get_first_nt_adjusted_coord(res: 'Residue') -> int:
            if len(res.exons) == 1 or res.codon[0].exon is res.codon[1].exon:
                return res.codon[0].coordinate
            else:  # if codon's primary exon is downstream, "abacus" the coord of the 1st nucleotide
                return res.codon[1].coordinate - 1


class GapResidue(Residue):
    def __init__(self, protein: 'Protein', position: int, upstream_exon: 'Exon', downstream_exon: 'Exon'):
        # if upstream_exon.transcript is not downstream_exon.transcript:
        #     raise ValueError(f'{upstream_exon} and {downstream_exon} belong to different transcripts')
        # if upstream_exon.transcript is not protein.orf.transcript:
        #     raise ValueError(f'{upstream_exon} and {protein} belong to different transcripts')
        super().__init__(protein, AminoAcid.GAP, position)
        self.upstream_exon = upstream_exon
        self.downstream_exon = downstream_exon
    
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
    def __init__(self, anchor: Optional['Residue'], other: Optional['Residue'], ttype: TranscriptLevelEvent):
        super().__init__(anchor, other)
        self.ttype = ttype
        self.ptype = None


class ProteinAlignment(Alignment):
    def __init__(self, anchor: 'Protein', other: 'Protein'):
        if anchor.orf.gene is not other.orf.gene:
            raise ValueError(f'{anchor} and {other} belong to different genes')
        if anchor.orf.transcript.strand is not other.orf.transcript.strand:
            raise ValueError(f'{anchor.orf.transcript} and {other.orf.transcript} are on different strands')
        super().__init__(anchor, other)
        
        anchor_stack = deque(anchor.residues)
        other_stack = deque(other.residues)
        anchor_res, other_res = None, None
        self.chain = []
        while anchor_stack or other_stack:
            if anchor_stack:
                anchor_res_coord = get_first_nt_adjusted_coord(anchor_stack[0])
            else:
                anchor_res_coord = inf
            if other_stack:
                other_res_coord = get_first_nt_adjusted_coord(other_stack[0])
            else:
                other_res_coord = inf
            if anchor.orf.transcript.strand is Strand.MINUS:
                if isfinite(anchor_res_coord):
                    anchor_res_coord = -anchor_res_coord
                if isfinite(other_res_coord):
                    other_res_coord = -other_res_coord

            if abs(anchor_res_coord - other_res_coord) < 2:
                anchor_res = anchor_stack.popleft()
                other_res = other_stack.popleft()
                res_align = ResidueAlignment(
                    anchor_res,
                    other_res,
                    TranscriptLevelEvent.MATCH if anchor_res.amino_acid is other_res.amino_acid else TranscriptLevelEvent.FRAMESHIFT
                )
            # TODO: set upstream and downstream exons for GapResidue
            elif anchor_res_coord > other_res_coord:
                gap_pos = anchor_res.position if anchor_res else 0
                other_res = other_stack.popleft()
                res_align = ResidueAlignment(
                    GapResidue(anchor, gap_pos, None, None),
                    other_res,
                    TranscriptLevelEvent.INSERTION
                )                
            elif anchor_res_coord < other_res_coord:
                gap_pos = other_res.position if other_res else 0
                anchor_res = anchor_stack.popleft()
                res_align = ResidueAlignment(
                    anchor_res,
                    GapResidue(other, gap_pos, None, None),
                    TranscriptLevelEvent.DELETION
                )
            assert res_align.anchor.protein is anchor
            assert res_align.other.protein is other
            self.chain.append(res_align)


    
