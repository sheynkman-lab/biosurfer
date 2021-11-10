from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from operator import attrgetter
from typing import TYPE_CHECKING, List, Optional, Tuple
from warnings import warn

from biosurfer.core.constants import Nucleobase, AminoAcid, Strand
from biosurfer.core.helpers import frozendataclass

if TYPE_CHECKING:
    from biosurfer.core.models.biomolecules import Chromosome, Gene, Transcript, Exon, ORF, Protein


# TODO: save memory usage with __slots__?
class Nucleotide:
    __slots__ = ('parent', 'coordinate', 'position', '_base', 'residue')

    def __init__(self, parent, coordinate: int, position: int) -> None:
        self.parent = parent
        self.coordinate = coordinate  # genomic coordinate
        self.position = position  # position within parent
        self._base = None
        self.residue = None  # associated Residue, if any
    
    def __repr__(self) -> str:
        return f'{self.parent.chromosome}:{self.coordinate}({self.parent.strand}){self.base}'
    
    @property
    def base(self) -> 'Nucleobase':
        if not self._base:
            self._base = Nucleobase(self.parent.sequence[self.position - 1])
        return self._base

    @property
    def chromosome(self) -> 'Chromosome':
        return self.parent.chromosome
    
    @property
    def strand(self) -> 'Strand':
        return self.parent.strand

    @property
    def gene(self) -> 'Gene':
        if isinstance(self.parent, Gene):
            return self.parent
        return self.parent.gene
    
    @property
    def exon(self) -> 'Exon':
        return self.parent.get_exon_containing_position(self.position)
    
    def __eq__(self, other: 'Nucleotide'):
        if not isinstance(other, Nucleotide):
            raise TypeError(f'Cannot compare Nucleotide with {type(other)}')
        return self.coordinate, self.position, self.base == other.coordinate, other.position, other.base


OptNucleotide = Optional['Nucleotide']
Codon = Tuple[OptNucleotide, OptNucleotide, OptNucleotide]

class Residue:
    __slots__ = ('protein', 'position', '_aa', 'codon')

    def __init__(self, protein: 'Protein', position: int) -> None:
        self._aa = None
        self.protein = protein
        self.position = position  # position within protein peptide sequence
        self.codon: Codon = (None, None, None)  # 3-tuple of associated Nucleotides; filled in later
    
    def __repr__(self) -> str:
        return f'{self.amino_acid}{self.position}'
    
    @property
    def amino_acid(self) -> 'AminoAcid':
        if not self._aa:
            self._aa = AminoAcid(self.protein.sequence[self.position - 1])
        return self._aa

    @property
    def codon_str(self) -> str:
        return ''.join(str(nt.base) for nt in self.codon)

    @property
    def exons(self) -> List['Exon']:
        return [nt.exon for nt in self.codon]
    
    @property
    def primary_exon(self) -> 'Exon':
        exons = Counter([nt.exon for nt in self.codon])
        return exons.most_common(1)[0][0]
    
    @property
    def junction(self) -> Optional['Junction']:
        exons = self.exons
        if len(exons) < 2:
            return None
        transcript = exons[0].transcript
        # this is a bit kludgy, may need to add properties to Exon class
        return transcript.junctions[exons[0].position - 1]
    
    @property
    def is_gap(self):
        return self.amino_acid is AminoAcid.GAP


@frozendataclass
class Junction:
    donor: int
    acceptor: int
    chromosome: 'Chromosome'
    strand: 'Strand'
    
    def __repr__(self) -> str:
        return f'{self.chromosome}({self.strand}):{self.donor}^{self.acceptor}'
    
    def __eq__(self, other: 'Junction') -> bool:
        if not isinstance(other, Junction):
            raise TypeError(f'Cannot compare Junction with {type(other)}')
        if self.chromosome is not other.chromosome:
            return False
        if self.strand is not other.strand:
            return False
        delta_donor = abs(self.donor - other.donor)
        delta_acceptor = abs(self.acceptor - other.acceptor)
        if delta_donor <= 2 and delta_acceptor <= 2 and (delta_donor != 0 or delta_acceptor != 0):
            warn(f'Possible off-by-one error for junctions {self} and {other}')
        if self.strand is Strand.MINUS and self.donor == other.acceptor and self.acceptor == other.donor:
            warn(f'Reversed coords for junctions {self} and {other}')
        return delta_donor == 0 and delta_acceptor == 0


class UTR(ABC):
    def __init__(self, orf: 'ORF', boundary_exon_index: int):
        self.orf = orf
        self.transcript: 'Transcript' = orf.transcript
        self._boundary_exon_index = boundary_exon_index
        self.transcript_start = None
        self.transcript_stop = None

    @property
    def length(self):
        return self.transcript_stop - self.transcript_start + 1

    @property
    def nucleotides(self) -> List['Nucleotide']:
        return self.transcript.nucleotides[self.transcript_start-1:self.transcript_stop]
    
    @property
    def sequence(self) -> str:
        return self.transcript.sequence[self.transcript_start-1:self.transcript_stop]
    
    @property
    def start(self) -> int:
        return self.transcript.nucleotides[self.transcript_start-1].coordinate
    
    @property
    def stop(self) -> int:
        return self.transcript.nucleotides[self.transcript_stop-1].coordinate
    
    @property
    @abstractmethod
    def exons(self):
        raise NotImplementedError

    @abstractmethod    
    def __repr__(self):
        raise NotImplementedError


class FivePrimeUTR(UTR):
    def __init__(self, orf: 'ORF', boundary_exon_index: int):
        super().__init__(orf, boundary_exon_index)
        self.transcript_start = 1
        self.transcript_stop = orf.transcript_start - 1

    @property
    def exons(self) -> List['Exon']:
        return self.transcript.exons[:self._boundary_exon_index+1]

    def __repr__(self):
        return f'{self.transcript}:utr5({self.transcript_start}-{self.transcript_stop})'


class ThreePrimeUTR(UTR):
    def __init__(self, orf: 'ORF', boundary_exon_index: int):
        super().__init__(orf, boundary_exon_index)
        self.transcript_start = orf.transcript_stop + 1
        self.transcript_stop = self.transcript.length
    
    @property
    def exons(self) -> List['Exon']:
        return self.transcript.exons[self._boundary_exon_index:]
    
    def __repr__(self):
        return f'{self.transcript}:utr3({self.transcript_start}-{self.transcript_stop})'
