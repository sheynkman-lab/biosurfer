from abc import ABC
from typing import Optional
from constants import ProteinLevelEvent, TranscriptLevelEvent
from models import Transcript, Exon, ORF, Nucleotide, Protein, Residue


class Alignment(ABC):
    def __init__(self, anchor, other):
        self.anchor = anchor
        self.other = other
    
    def __repr__(self):
        return f'{self.anchor}|{self.other}'


class ResidueAlignment(Alignment):
    def __init__(self, anchor: Optional['Residue'], other: Optional['Residue'], ptype: ProteinLevelEvent, ttype: TranscriptLevelEvent):
        super().__init__(anchor, other)


class ProteinAlignment(Alignment):
    def __init__(self, anchor: 'Protein', other: 'Protein'):
        super().__init__(anchor, other)
        
        def get_abacused_coords():
            pass

    
