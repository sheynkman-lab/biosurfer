from abc import ABC, abstractmethod
from collections import Counter
from functools import cached_property
from operator import attrgetter
from typing import Dict, Iterable, List, Optional, Tuple, Type
from warnings import warn

from Bio.Seq import Seq
from biosurfer.core.constants import (APPRIS, SQANTI, AminoAcid, Nucleobase,
                                      Strand)
from biosurfer.core.helpers import BisectDict, frozendataclass
from sqlalchemy import (Boolean, Column, Enum, ForeignKey, Integer, String,
                        create_engine)
from sqlalchemy.ext.declarative import (DeclarativeMeta, declarative_base,
                                        has_inherited_table)
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import reconstructor, relationship
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql import func, select


class BaseMeta(DeclarativeMeta):
    _session = None

    def __init__(cls, classname, bases, dict_, **kwargs):
        cls.__tablename__ = cls.__name__.lower() if not has_inherited_table(cls) else None
        DeclarativeMeta.__init__(cls, classname, bases, dict_, **kwargs)

    @property
    def session(cls):
        if BaseMeta._session is None:
            raise AttributeError('Base.session has not been set')
        return BaseMeta._session
    
    @session.setter
    def session(cls, session):
        BaseMeta._session = session


Base = declarative_base(metaclass=BaseMeta)


class NameMixin:
    name = Column(String, index=True)

    @classmethod
    def from_name(cls: Type['Base'], name: str, unique: bool = True):
        statement = select(cls).where(cls.name == name)
        result = cls.session.execute(statement).scalars()
        if unique:
            try:
                return result.one()
            except NoResultFound:
                return None
        else:
            return result.all()

    @classmethod
    def from_names(cls: Type['Base'], names: Iterable[str]):
        statement = select(cls).where(cls.name.in_(names))
        return {inst.name: inst for inst in cls.session.execute(statement).scalars()}


class AccessionMixin:
    accession = Column(String, primary_key=True, index=True)

    @classmethod
    def from_accession(cls: Type['Base'], accession: str):
        statement = select(cls).where(cls.accession == accession)
        result = cls.session.execute(statement).scalars()
        try:
            return result.one()
        except NoResultFound:
            return None
    
    @classmethod
    def from_accessions(cls: Type['Base'], accessions: Iterable[str]):
        statement = select(cls).where(cls.accession.in_(accessions))
        return {inst.name: inst for inst in cls.session.execute(statement).scalars()}


class Chromosome(Base, NameMixin):
    name = Column(String, primary_key=True)
    genes = relationship('Gene', back_populates='chromosome')

    def __repr__(self) -> str:
        return self.name


class Gene(Base, NameMixin, AccessionMixin):
    strand = Column(Enum(Strand))
    chromosome_id = Column(String, ForeignKey('chromosome.name'))
    chromosome = relationship('Chromosome', back_populates='genes')
    transcripts = relationship(
        'Transcript',
        back_populates='gene',
        order_by='Transcript.name',
        lazy='selectin'  # always load transcripts along with gene
    )

    def __repr__(self) -> str:
        return self.name

    # FIXME: make this work in SQL queries
    @property
    def start(self) -> int:
        return min(exon.start for transcript in self.transcripts for exon in transcript.exons)
    
    @property
    def stop(self) -> int:
        return max(exon.stop for transcript in self.transcripts for exon in transcript.exons)
    


class Transcript(Base, NameMixin, AccessionMixin):
    strand = Column(Enum(Strand))
    type = Column(String)
    sequence = Column(String)
    gene_id = Column(String, ForeignKey('gene.accession'))
    gene = relationship('Gene', back_populates='transcripts')
    exons = relationship(
        'Exon',
        order_by='Exon.transcript_start',
        collection_class=ordering_list('position', count_from=1),
        back_populates='transcript',
        uselist=True,
        lazy='selectin'  # always load exons along with transcript
    )
    orfs = relationship(
        'ORF',
        order_by='ORF.transcript_start',
        back_populates='transcript',
        uselist=True
    )

    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': 'transcript'
    }

    # def __init__(self, **kwargs):
    #     super().__init__(**kwargs)
    #     self.init_on_load()

    # @reconstructor
    # def init_on_load(self):
    #     self._nucleotide_mapping: Dict[int, 'Nucleotide'] = dict()

    # The reason we use cached properties here instead of setting things up in __init__ or init_on_load
    # is to make sure the ORM eagerly loads all exons first.
    @cached_property
    def _exon_mapping(self) -> BisectDict:
        return BisectDict((exon.transcript_stop+1, i) for i, exon in enumerate(self.exons))

    @cached_property
    def _junction_mapping(self) -> Dict['Junction', Tuple['Exon', 'Exon']]:
        mapping = dict()
        for i in range(1, len(self.exons)):
            up_exon = self.exons[i-1]
            down_exon = self.exons[i]

            donor = up_exon.nucleotides[-1].coordinate
            acceptor = down_exon.nucleotides[0].coordinate

            junction = Junction(donor, acceptor, self.chromosome, self.strand)
            mapping[junction] = (up_exon, down_exon)
        return mapping
    
    @cached_property
    def nucleotides(self):
        if not self.sequence:
            raise AttributeError(f'{self.name} has no sequence')
        # assert sum(exon.length for exon in self.exons) == self.length
        nucleotides = []
        self._nucleotide_mapping: Dict[int, 'Nucleotide'] = dict()
        if self.exons:
            i = 0
            for exon in self.exons:
                coords = range(exon.start, exon.stop+1)
                if self.strand is Strand.MINUS:
                    coords = reversed(coords)
                for coord in coords:
                    nt = Nucleotide(self, coord, i+1, self.sequence[i])
                    nucleotides.append(nt)
                    self._nucleotide_mapping[coord] = nt
                    i += 1
        else:
            nucleotides = [Nucleotide(self, None, i, base) for i, base in enumerate(self.sequence, start=1)]
        return nucleotides

    def __repr__(self) -> str:
        return self.name

    # FIXME: make this work in SQL queries
    @property
    def start(self) -> int:
        return min(exon.start for exon in self.exons)
    
    @property
    def stop(self) -> int:
        return max(exon.stop for exon in self.exons)
    
    @hybrid_property
    def length(self):
        return len(self.sequence)
    
    @length.expression
    def length(cls):
        return func.length(cls.sequence)

    @property
    def chromosome(self) -> 'Chromosome':
        return self.gene.chromosome
    
    @property
    def primary_orf(self) -> Optional['ORF']:
        if not self.orfs:
            return None
        return max(self.orfs, key=attrgetter('length'))

    @property
    def protein(self) -> Optional['Protein']:
        """Get the "primary" protein produced by this transcript, if it exists."""
        return self.primary_orf.protein if self.primary_orf else None
    
    @property
    def junctions(self):
        return list(self._junction_mapping.keys())

    # These methods may seem redundant, but the idea is to keep the publicly accessible interface separate from the implementation details
    def get_exon_containing_position(self, position: int) -> 'Exon':
        """Given a position (1-based) within the transcript's nucleotide sequence, return the exon containing that position."""
        return self.exons[self._exon_mapping[position]]
    
    def get_exon_index_containing_position(self, position: int) -> int:
        """Given a position (1-based) within the transcript's nucleotide sequence, return the index (0-based) of the exon containing that position."""
        return self._exon_mapping[position]
    
    def get_nucleotide_from_coordinate(self, coordinate: int) -> 'Nucleotide':
        """Given a genomic coordinate (1-based) included in the transcript, return the Nucleotide object corresponding to that coordinate."""
        if coordinate in self._nucleotide_mapping:
            return self._nucleotide_mapping[coordinate]
        else:
            return None

    def contains_coordinate(self, coordinate: int) -> bool:
        """Given a genomic coordinate (1-based), return whether or not the transcript contains the nucleotide at that coordinate."""
        return coordinate in self._nucleotide_mapping
    
    def get_exons_from_junction(self, junction: 'Junction') -> Tuple['Exon', 'Exon']:
        try:
            return self._junction_mapping[junction]
        except KeyError as e:
            raise KeyError(f'{self} does not use junction {junction}') from e


class GencodeTranscript(Transcript):
    appris = Column(Enum(APPRIS))
    start_nf = Column(Boolean)
    end_nf = Column(Boolean)
    pacbio = relationship(
        'PacBioTranscript',
        back_populates = 'gencode',
        uselist = True
    )
    
    def __init__(self, **kwargs):
        if 'strand' in kwargs:
            kwargs['strand'] = Strand.from_symbol(kwargs['strand'])
        super().__init__(**kwargs)
    
    __mapper_args__ = {
        'polymorphic_identity': 'gencodetranscript'
    }

    @hybrid_property
    def basic(self):
        return ~(self.start_nf | self.end_nf)


class PacBioTranscript(Transcript):
    sqanti = Column(Enum(SQANTI))
    gencode_id = Column(String, ForeignKey('transcript.accession'))
    gencode = relationship(
        'GencodeTranscript',
        back_populates = 'pacbio',
        uselist = False,
        remote_side = [Transcript.accession]
    )

    __mapper_args__ = {
        'polymorphic_identity': 'pacbiotranscript'
    }


class Exon(Base, AccessionMixin):
    type = Column(String)
    position = Column(Integer)  # exon ordinal within parent transcript
    # genomic coordinates
    start = Column(Integer)
    stop = Column(Integer)
    # transcript coordinates
    transcript_start = Column(Integer)
    transcript_stop = Column(Integer)
    # sequence = Column(String, default='')
    transcript_id = Column(String, ForeignKey('transcript.accession'), primary_key=True)
    transcript = relationship(
        'Transcript', 
        back_populates='exons'
    )
    
    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': 'exon'
    }

    def __repr__(self) -> str:
        return f'{self.transcript}:exon{self.position}'
  
    @hybrid_property
    def length(self):
        return self.stop - self.start + 1
    
    @property
    def gene(self):
        return self.transcript.gene
        
    @property
    def chromosome(self):
        return self.gene.chromosome
    
    @property
    def strand(self) -> 'Strand':
        return self.transcript.strand
    
    @property
    def sequence(self):
        return self.transcript.sequence[self.transcript_start-1:self.transcript_stop]

    @property
    def nucleotides(self):
        return self.transcript.nucleotides[self.transcript_start-1:self.transcript_stop]

    @property
    def coding_nucleotides(self):
        return [nt for nt in self.nucleotides if nt.residue]


class GencodeExon(Exon):
    __mapper_args__ = {
        'polymorphic_identity': 'gencodeexon'
    }


class PacBioExon(Exon):
    __mapper_args__ = {
        'polymorphic_identity': 'pacbioexon'
    }


# TODO: save memory usage with __slots__?
class Nucleotide:
    def __init__(self, parent, coordinate: int, position: int, base: str) -> None:
        self.parent = parent
        self.coordinate = coordinate  # genomic coordinate
        self.position = position  # position within parent
        self.base = Nucleobase(base)
        self.residue = None  # associated Residue, if any
    
    def __repr__(self) -> str:
        return f'{self.parent.chromosome}:{self.coordinate}({self.parent.strand}){self.base}'
    
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


class ORF(Base):
    # genomic coordinates
    # start = Column(Integer)  
    # stop = Column(Integer)  
    # transcript coordinates
    transcript_start = Column(Integer, primary_key=True)
    transcript_stop = Column(Integer, primary_key=True)
    transcript_id = Column(String, ForeignKey('transcript.accession'), primary_key=True)
    transcript = relationship(
        'Transcript', 
        back_populates='orfs'
    )
    protein_id = Column(String, ForeignKey('protein.accession'))
    protein = relationship(
        'Protein',
        back_populates='orf',
        uselist=False
    )

    # __table_args__ = (UniqueConstraint(transcript_id, transcript_start, transcript_stop, name='_transcript_orf_pair'),)

    @reconstructor
    def init_on_load(self):
        self._first_exon_index = self.transcript.get_exon_index_containing_position(self.transcript_start)
        try:
            self._last_exon_index = self.transcript.get_exon_index_containing_position(self.transcript_stop)
        except KeyError as e:
            warn(f'KeyError: {e} when initializing {self}')
            self._last_exon_index = len(self.transcript.exons) - 1

        utr5_boundary_exon_index = self._first_exon_index
        utr3_boundary_exon_index = self._last_exon_index
        if self.transcript.exons[utr5_boundary_exon_index].transcript_start == self.transcript_start:
            utr5_boundary_exon_index -= 1
        if self.transcript.exons[utr3_boundary_exon_index].transcript_stop == self.transcript_stop:
            utr3_boundary_exon_index += 1

        if utr5_boundary_exon_index >= 0:
            self.utr5 = FivePrimeUTR(self, utr5_boundary_exon_index)
        else:
            self.utr5 = None
        if utr3_boundary_exon_index < len(self.transcript.exons):
            self.utr3 = ThreePrimeUTR(self, utr3_boundary_exon_index)
        else:
            self.utr3 = None

        # ORFs with stop codons at least 50 bp upstream of the last splice site in the mature transcript
        # (i.e. the beginning of the last exon) are considered candidates for nonsense-mediated decay (NMD)
        last_junction = self.transcript.exons[-1].transcript_start
        self.nmd = last_junction - self.transcript_stop >= 50

    def __repr__(self) -> str:
        return f'{self.transcript}:orf({self.transcript_start}-{self.transcript_stop})'

    # FIXME: make this work in SQL queries
    @property
    def start(self):
        if self.transcript.strand is Strand.PLUS:
            transcript_coord = self.transcript_start
        elif self.transcript.strand is Strand.MINUS:
            transcript_coord = self.transcript_stop
        else:
            return None
        return self.transcript.nucleotides[transcript_coord - 1].coordinate

    @property
    def stop(self):
        if self.transcript.strand is Strand.PLUS:
            transcript_coord = self.transcript_stop
        elif self.transcript.strand is Strand.MINUS:
            transcript_coord = self.transcript_start
        else:
            return None
        return self.transcript.nucleotides[transcript_coord - 1].coordinate
    
    @hybrid_property
    def length(self) -> int:
        return self.transcript_stop - self.transcript_start + 1

    @property
    def sequence(self) -> str:
        return self.transcript.sequence[self.transcript_start - 1:self.transcript_stop]
    
    @property
    def nucleotides(self) -> List['Nucleotide']:
        return self.transcript.nucleotides[self.transcript_start - 1:self.transcript_stop]
    
    @property
    def gene(self) -> 'Gene':
        return self.transcript.gene
    
    @property
    def exons(self) -> List['Exon']:
        return self.transcript.exons[self._first_exon_index:self._last_exon_index+1]
    
    @property
    def junctions(self) -> List['Junction']:
        return self.transcript.junctions[self._first_exon_index:self._last_exon_index]
    
    def _link_aa_to_nt(self, residue_list):
        aa_sequence = Seq(self.protein.sequence)
        nt_sequence = Seq(self.sequence)
        translation = nt_sequence.translate(to_stop=True)
        aa_match_index = aa_sequence.find(translation)
        if aa_match_index == -1:
            warn(
                f'Could not match amino acid sequence to nucleotide sequence of {self}'
            )
            return
        
        nt_match_index = aa_match_index*3
        nt_list = self.nucleotides[nt_match_index:]
        for i, aa in enumerate(residue_list):
            aa.codon = tuple(nt_list[3*i:3*i + 3])
            for nt in aa.codon:
                nt.residue = aa


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


class Protein(Base, AccessionMixin):
    sequence = Column(String)
    orf = relationship(
        'ORF',
        back_populates = 'protein',
        uselist = False
    )
    features = relationship(
        'ProteinFeature',
        back_populates = 'protein',
        uselist = True
    )

    # def __init__(self, **kwargs):
    #     super().__init__(**kwargs)
    #     self.residues = []

    # @reconstructor
    # def init_on_load(self):
    #     self.residues = [Residue(self, aa, i) for i, aa in enumerate(self.sequence + '*', start=1)]
    #     if self.orf and self.orf.nucleotides:
    #         self._link_aa_to_orf_nt()
    
    @cached_property
    def residues(self):
        if not self.sequence.endswith('*'):
            self.sequence += '*'
        _residues = [Residue(self, aa, i) for i, aa in enumerate(self.sequence, start=1)]
        if self.orf:
            self.orf._link_aa_to_nt(_residues)
        return _residues

    def __repr__(self):
        return f'{self.orf.transcript}:protein'
    
    @property
    def gene(self):
        return self.orf.transcript.gene
    
    @property
    def transcript(self):
        return self.orf.transcript
    
    @hybrid_property
    def length(self):
        return len(self.sequence)


class ProteinFeature(Base):
    type = Column(String)
    name = Column(String)
    accession = Column(String, primary_key=True)
    protein_id = Column(String, ForeignKey('protein.accession'), primary_key=True)
    protein_start = Column(Integer, primary_key=True)
    protein_stop = Column(Integer, primary_key=True)
    protein = relationship(
        'Protein',
        back_populates = 'features',
        uselist = False
    )

    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': 'feature'
    }
    
    def __repr__(self):
        return f'{self.protein}:{self.name}({self.protein_start}-{self.protein_stop})'

    @hybrid_property
    def length(self):
        return self.protein_stop - self.protein_start + 1

    @property
    def sequence(self) -> str:
        return self.protein.sequence[self.protein_start-1:self.protein_stop]
    
    @property
    def residues(self) -> List['Residue']:
        return self.protein.residues[self.protein_start-1:self.protein_stop]


OptNucleotide = Optional['Nucleotide']
Codon = Tuple[OptNucleotide, OptNucleotide, OptNucleotide]

class Residue:
    def __init__(self, protein: 'Protein', amino_acid: str, position: int) -> None:
        self.amino_acid = AminoAcid(amino_acid)
        self.protein = protein
        self.position = position  # position within protein peptide sequence
        self.codon: Codon = (None, None, None)  # 3-tuple of associated Nucleotides; filled in later
    
    def __repr__(self) -> str:
        return f'{self.amino_acid}{self.position}'
    
    @property
    def codon_str(self) -> str:
        return ''.join(str(nt.base) for nt in self.codon)

    @property
    def exons(self) -> List['Exon']:
        # TODO: is sorting necessary here?
        return sorted({nt.exon for nt in self.codon}, key=attrgetter('position'))
    
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
