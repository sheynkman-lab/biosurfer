from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from operator import attrgetter
from typing import TYPE_CHECKING, List, Optional, Tuple, Union
from warnings import warn

from attrs import define, frozen, field, evolve, validators
from biosurfer.core.constants import Nucleobase, AminoAcid, Strand

if TYPE_CHECKING:
    from biosurfer.core.models.biomolecules import Chromosome, Gene, Transcript, Exon, ORF, Protein


@frozen(hash=True)
class Position:
    chromosome: str  # TODO: convert to dynamic Enum?
    strand: 'Strand' = field(validator=validators.instance_of(Strand))
    coordinate: int = field(validator=validators.gt(0))

    def __repr__(self):
        return f'{self.chromosome}({self.strand}):{self.coordinate}'

    def _is_comparable(self, other: 'Position'):
        return isinstance(other, Position) and (self.chromosome, self.strand) == (other.chromosome, other.strand)

    def __lt__(self, other: 'Position'):
        if not self._is_comparable(other):
            return NotImplemented
        elif self.strand is Strand.MINUS:
            return self.coordinate > other.coordinate
        else:
            return self.coordinate < other.coordinate
    
    def __le__(self, other: 'Position'):
        if not self._is_comparable(other):
            return NotImplemented
        elif self.strand is Strand.MINUS:
            return self.coordinate >= other.coordinate
        else:
            return self.coordinate <= other.coordinate
    
    def __gt__(self, other: 'Position'):
        if not self._is_comparable(other):
            return NotImplemented
        elif self.strand is Strand.MINUS:
            return self.coordinate < other.coordinate
        else:
            return self.coordinate > other.coordinate
    
    def __ge__(self, other: 'Position'):
        if not self._is_comparable(other):
            return NotImplemented
        elif self.strand is Strand.MINUS:
            return self.coordinate <= other.coordinate
        else:
            return self.coordinate >= other.coordinate
    
    def __add__(self, offset: int):
        if not isinstance(offset, int):
            return NotImplemented
        else:
            if self.strand is Strand.MINUS:
                offset = -offset
            return evolve(self, coordinate=self.coordinate + offset)
    
    def __radd__(self, offset: int):
        return self.__add__(offset)

    def __sub__(self, other: Union['Position', int]):
        if isinstance(other, int):
            return self.__add__(-other)
        elif not self._is_comparable(other):
            return NotImplemented
        delta = self.coordinate - other.coordinate
        if self.strand is Strand.MINUS:
            delta = -delta
        return delta


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
        return list(set(nt.exon for nt in self.codon))
    
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


@frozen(hash=True)
class GenomeRange:
    begin: 'Position' = field(validator=validators.instance_of(Position))
    end: 'Position' = field(validator=validators.instance_of(Position))
    
    def __attrs_post_init__(self):
        if self.begin > self.end:
            raise ValueError(f'{self.begin} is downstream of {self.end}')

    @property
    def chromosome(self):
        return self.begin.chromosome
    
    @property
    def strand(self):
        return self.begin.strand

    def __repr__(self):
        return f'{self.begin.chromosome}({self.begin.strand}):{self.begin.coordinate}^{self.end.coordinate}'
    
    def __eq__(self, other: 'GenomeRange'):
        if not isinstance(other, GenomeRange):
            return NotImplemented
        delta_begin = abs(self.begin - other.begin)
        delta_end = abs(self.end - other.end)
        if delta_begin <= 2 and delta_end <= 2 and (delta_begin != 0 or delta_end != 0):
            warn(f'Possible off-by-one error for ranges {self} and {other}')
        return delta_begin == 0 and delta_end == 0
    
    def __and__(self, other: 'GenomeRange'):
        if not isinstance(other, GenomeRange):
            return NotImplemented
        begin = max(self.begin, other.begin)
        end = min(self.end, other.end)
        return evolve(self, begin=begin, end=end) if begin <= end else None

    def __or__(self, other: 'GenomeRange'):
        if not isinstance(other, GenomeRange):
            return NotImplemented
        begin = min(self.begin, other.begin)
        end = max(self.end, other.end)
        return evolve(self, begin=begin, end=end) if begin <= end else None

    @property
    def length(self) -> int:
        return (self.end - self.begin) + 1

    def as_tuple(self):
        return self.chromosome, self.strand, self.begin.coordinate, self.end.coordinate

    @classmethod
    def from_coordinates(cls, chromosome: str, strand: 'Strand', begin: int, end: int):
        return cls(Position(chromosome, strand, begin), Position(chromosome, strand, end))


@frozen(hash=True)
class Junction:
    range: 'GenomeRange' = field(validator=validators.instance_of(GenomeRange))

    @property
    def donor(self):
        return self.range.begin
    
    @property
    def acceptor(self):
        return self.range.end

    def __repr__(self):
        return f'{self.range.chromosome}({self.range.strand}):{self.donor.coordinate}^{self.acceptor.coordinate}'
    
    def __eq__(self, other: 'Junction'):
        if not isinstance(other, Junction):
            return NotImplemented
        return self.range == other.range
    
    def __and__(self, other: 'Junction'):
        if not isinstance(other, Junction):
            return NotImplemented
        intersection = self.range & other.range
        return Junction(intersection) if intersection else None

    def __or__(self, other: 'Junction'):
        if not isinstance(other, Junction):
            return NotImplemented
        union = self.range | other.range
        return Junction(union) if union else None

    @property
    def length(self):
        return self.range.length

    def as_tuple(self):
        return self.range.as_tuple()

    @classmethod
    def from_coordinates(cls, chromosome: str, strand: 'Strand', donor: int, acceptor: int):
        return Junction(GenomeRange.from_coordinates(chromosome, strand, donor, acceptor))

    @classmethod
    def from_splice_sites(cls, donor: 'Position', acceptor: 'Position'):
        return Junction(GenomeRange(donor, acceptor))


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
