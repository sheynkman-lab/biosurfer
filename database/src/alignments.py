from abc import ABC
from collections import deque
from math import inf, isfinite
from typing import Optional

from constants import AminoAcid, ProteinLevelEvent, TranscriptLevelEvent
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
        anchor_gap, other_gap = None, None
        self.chain = []
        while anchor_stack or other_stack:
            if anchor_stack:
                # anchor_res_coord = get_first_nt_adjusted_coord(anchor_stack[0])
                anchor_res_coord = anchor_stack[0].codon[0].coordinate
            else:
                anchor_res_coord = inf
            if other_stack:
                # other_res_coord = get_first_nt_adjusted_coord(other_stack[0])
                other_res_coord = other_stack[0].codon[0].coordinate
            else:
                other_res_coord = inf
            if anchor.orf.transcript.strand is Strand.MINUS:
                if isfinite(anchor_res_coord):
                    anchor_res_coord = -anchor_res_coord
                if isfinite(other_res_coord):
                    other_res_coord = -other_res_coord
            
            coord_diff = anchor_res_coord - other_res_coord
            if abs(coord_diff) < 2:
                anchor_gap, other_gap = None, None
                anchor_res = anchor_stack.popleft()
                other_res = other_stack.popleft()
                if coord_diff != 0:
                    event_type = TranscriptLevelEvent.FRAMESHIFT
                elif anchor_res.amino_acid is not other_res.amino_acid:
                    event_type = TranscriptLevelEvent.SPLIT
                else:
                    event_type = TranscriptLevelEvent.MATCH
                res_align = ResidueAlignment(anchor_res, other_res, event_type)
            # TODO: set upstream and downstream exons for GapResidue
            elif coord_diff > 0:
                other_gap = None
                if not anchor_gap:
                    gap_pos = anchor_res.position if anchor_res else 0
                    anchor_gap = GapResidue(anchor, gap_pos, None, None)
                other_res = other_stack.popleft()
                res_align = ResidueAlignment(anchor_gap, other_res, TranscriptLevelEvent.INSERTION)                
            elif coord_diff < 0:
                anchor_gap = None
                if not other_gap:
                    gap_pos = other_res.position if other_res else 0
                    other_gap = GapResidue(other, gap_pos, None, None)
                anchor_res = anchor_stack.popleft()
                res_align = ResidueAlignment(anchor_res, other_gap, TranscriptLevelEvent.DELETION)
            assert res_align.anchor.protein is anchor
            assert res_align.other.protein is other
            self.chain.append(res_align)
    
    @property
    def full(self):
        anchor_str = ''.join(str(res.anchor.amino_acid) for res in self.chain)
        other_str = ''.join(str(res.other.amino_acid) for res in self.chain)
        ttype_str = ''.join(str(res.ttype) for res in self.chain)
        return anchor_str + '\n' + ttype_str + '\n' + other_str
    


    
