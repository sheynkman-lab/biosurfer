from collections import Counter
from operator import attrgetter
from typing import List, Optional, Tuple
from warnings import warn

from Bio.Seq import Seq
from inscripta.biocantor.location.location_impl import (CompoundInterval,
                                                        SingleInterval)
from sqlalchemy import Boolean, Column, Enum, ForeignKey, Integer, String, select
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import reconstructor, relationship
from sqlalchemy.orm.exc import NoResultFound

from constants import APPRIS, AminoAcid, Nucleobase, Strand, UTRType
from database import Base, db_session
from helpers import BisectDict


class Chromosome(Base):
    __tablename__ = 'chromosome'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    genes = relationship('Gene', back_populates='chromosome')

    def __repr__(self) -> str:
        return self.name


class Gene(Base):
    __tablename__ = 'gene'

    id = Column(Integer, primary_key=True)
    accession = Column(String)
    name = Column(String)
    chromosome_id = Column(Integer,ForeignKey('chromosome.id'))
    transcripts = relationship('Transcript', back_populates='gene', order_by='Transcript.name')
    chromosome = relationship('Chromosome', back_populates='genes')
    def __repr__(self) -> str:
        return self.name
    
    @classmethod
    def from_name(cls, name: str):
        statement = select(Gene).filter(Gene.name == name)
        try:
            return db_session.execute(statement).one()[Gene]
        except NoResultFound:
            return None
    
    @hybrid_property
    def strand(self) -> Strand:
        strands = Counter(transcript.strand for transcript in self.transcripts)
        return strands.most_common(1)[0][0]


class Transcript(Base):
    __tablename__ = 'transcript'
    id = Column(Integer, primary_key=True)
    accession = Column(String)
    name = Column(String)
    strand = Column(Enum(Strand))
    type = Column(String)
    sequence = Column(String)
    gene_id = Column(Integer, ForeignKey('gene.id'))
    gene = relationship('Gene', back_populates='transcripts')
    exons = relationship(
        'Exon',
        order_by='Exon.position',
        collection_class=ordering_list('position', count_from=1),
        back_populates='transcript',
        uselist=True
    )
    orfs = relationship(
        'ORF',
        order_by='ORF.start',
        back_populates='transcript',
        uselist=True)
    # TODO: add start and stop attrs

    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': 'transcript'
    }

    def _init_inner(self):
        self._nucleotides = None
        self._nucleotide_mapping = None

    def __init__(self):
        self._init_inner()

    @reconstructor
    def init_on_load(self):
        self._init_inner()
        self._exon_mapping = BisectDict((exon.transcript_stop+1, i) for i, exon in enumerate(self.exons))
    
    @hybrid_property
    def nucleotides(self):
        if not self.sequence:
            raise AttributeError(f'{self.name} has no sequence')
        elif self._nucleotides is None:
            assert sum(exon.stop - exon.start + 1 for exon in self.exons) == len(self.sequence)
            self._nucleotides = []
            self._nucleotide_mapping = dict()
            i = 0
            for exon in self.exons:
                coords = range(exon.start, exon.stop+1)
                if self.strand is Strand.MINUS:
                    coords = reversed(coords)
                for coord in coords:
                    nt = Nucleotide(self, coord, i+1, self.sequence[i])
                    self._nucleotides.append(nt)
                    self._nucleotide_mapping[coord] = nt
                    i += 1
        return self._nucleotides

    
    def __repr__(self) -> str:
        return self.name

    @hybrid_property
    def _location(self):
        exon = self.exons[0]
        loc = exon._location
        for exon in self.exons[1:]:
            loc = loc.union(exon._location)
        return loc

    @hybrid_property
    def start(self):
        return min(exon.start for exon in self.exons)
    
    @hybrid_property
    def stop(self):
        return max(exon.stop for exon in self.exons)
    
    @hybrid_property
    def length(self):
        return len(self.sequence)

    @hybrid_property
    def chromosome(self):
        return self.gene.chromosome
    
    @hybrid_property
    def protein(self):
        """Get the "primary" protein produced by this transcript, if it exists."""
        # TODO: implement this
        raise NotImplementedError
    
    # These methods may seem redundant, but the idea is to keep the publicly accessible interface separate from the implementation details
    def get_exon_containing_position(self, position: int) -> 'Exon':
        """Given a position (1-based) within the transcript's nucleotide sequence, return the exon containing that position."""
        return self.exons[self._exon_mapping[position]]
    
    def get_exon_index_containing_position(self, position: int) -> int:
        """Given a position (1-based) within the transcript's nucleotide sequence, return the number of the exon containing that position."""
        return self._exon_mapping[position]
    
    def get_nucleotide_from_coordinate(self, coordinate: int) -> 'Nucleotide':
        """Given a genomic coordinate (1-based) included in the transcript, return the Nucleotide object corresponding to that coordinate."""
        if coordinate in self._nucleotide_mapping:
            return self._nucleotide_mapping[coordinate]
        else:
            return None
    
    @classmethod
    def from_name(cls, name: str):
        statement = select(Transcript).filter(Transcript.name == name)
        try:
            return db_session.execute(statement).one()[Transcript]
        except NoResultFound:
            return None


class GencodeTranscript(Transcript):
    appris = Column(Enum(APPRIS))
    start_nf = Column(Boolean)
    end_nf = Column(Boolean)
    
    def __init__(self, *, accession, name, strand, appris, start_nf, end_nf):
        super().__init__()
        self.accession = accession
        self.name = name
        self.strand = Strand.from_symbol(strand)
        self.appris = appris
        self.start_nf = start_nf
        self.end_nf = end_nf
    
    __mapper_args__ = {
        'polymorphic_identity': 'gencode_transcript'
    }

    @hybrid_property
    def basic(self):
        return not (self.start_nf or self.end_nf)


class Exon(Base):
    __tablename__ = 'exon'
    id = Column(Integer, primary_key=True)
    accession = Column(String)
    type = Column(String)
    position = Column(Integer)  # exon ordinal within parent transcript
    # genomic coordinates
    start = Column(Integer)
    stop = Column(Integer)
    # transcript coordinates
    transcript_start = Column(Integer)
    transcript_stop = Column(Integer)
    # sequence = Column(String, default='')
    transcript_id = Column(Integer, ForeignKey('transcript.id'))
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
    def _location(self):
        return SingleInterval(self.transcript_start-1, self.transcript_stop, Strand.PLUS)
    
    @hybrid_property
    def length(self):
        return self.stop - self.start + 1
    
    @hybrid_property
    def gene(self):
        return self.transcript.gene
        
    @hybrid_property
    def chromosome(self):
        return self.gene.chromosome
    
    @hybrid_property
    def strand(self):
        return self.transcript.strand
    
    @hybrid_property
    def sequence(self):
        return self.transcript.sequence[self.transcript_start-1:self.transcript_stop]

    @hybrid_property
    def nucleotides(self):
        return self.transcript.nucleotides[self.transcript_start-1:self.transcript_stop]

    @hybrid_property
    def coding_nucleotides(self):
        return [nt for nt in self.nucleotides if nt.residue]


class GencodeExon(Exon):
    __mapper_args__ = {
        'polymorphic_identity': 'gencode_exon'
    }

    def __init__(self, *, accession, start, stop, transcript):
        self.accession = accession
        self.start = start
        self.stop = stop
        self.transcript = transcript


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
    def gene(self):
        if isinstance(self.parent, Gene):
            return self.parent
        return self.parent.gene
    
    @property
    def exon(self) -> 'Exon':
        return self.parent.get_exon_containing_position(self.position)
    

class ORF(Base):
    __tablename__ = 'orf'
    id = Column(Integer, primary_key=True)
    # genomic coordinates
    # TODO: pull these from first and last exon
    start = Column(Integer)  
    stop = Column(Integer)  
    # transcript coordinates
    transcript_start = Column(Integer)
    transcript_stop = Column(Integer)

    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    transcript = relationship(
        'Transcript', 
        back_populates='orfs'
    )
    protein = relationship(
        'Protein',
        back_populates='orf',
        uselist=False
    )

    @reconstructor
    def init_on_load(self):
        self._first_exon_index = self.transcript.get_exon_index_containing_position(self.transcript_start)
        self._last_exon_index = self.transcript.get_exon_index_containing_position(self.transcript_stop)
        self.utr5 = UTR(self, UTRType.FIVE_PRIME, self._first_exon_index)
        self.utr3 = UTR(self, UTRType.THREE_PRIME, self._last_exon_index)

    def __repr__(self) -> str:
        return f'{self.transcript}:orf({self.transcript_start}-{self.transcript_stop})'

    @hybrid_property
    def sequence(self):
        return self.transcript.sequence[self.transcript_start - 1:self.transcript_stop]
    
    @hybrid_property
    def nucleotides(self):
        return self.transcript.nucleotides[self.transcript_start - 1:self.transcript_stop]
    
    @hybrid_property
    def gene(self):
        return self.transcript.gene
    
    @hybrid_property
    def exons(self):
        return self.transcript.exons[self._first_exon_index:self._last_exon_index+1]

    @hybrid_property
    def _location(self):
        return SingleInterval(self.transcript_start-1, self.transcript_stop, Strand.PLUS)


class UTR:
    def __init__(self, orf: 'ORF', type: 'UTRType', boundary_exon_index: int):
        self.orf = orf
        self.transcript: 'Transcript' = orf.transcript
        self.type = type
        self._boundary_exon_index = boundary_exon_index

        if type is UTRType.FIVE_PRIME:
            self.transcript_start = 1
            self.transcript_stop = orf.transcript_start - 1
        elif type is UTRType.THREE_PRIME:
            self.transcript_start = orf.transcript_stop + 1
            self.transcript_stop = self.transcript.length
        self.length = self.transcript_stop - self.transcript_start + 1

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
    def exons(self) -> List['Exon']:
        if self.type is UTRType.FIVE_PRIME:
            return self.transcript.exons[:self._boundary_exon_index+1]
        elif self.type is UTRType.THREE_PRIME:
            return self.transcript.exons[self._boundary_exon_index:]
    
    def __repr__(self):
        return f'{self.transcript}:{self.type}({self.transcript_start}-{self.transcript_stop})'


class Protein(Base):
    __tablename__ = 'protein'
    id = Column(Integer, primary_key=True)
    
    sequence = Column(String)
    
    orf_id = Column(Integer, ForeignKey('orf.id'))
    orf = relationship(
        'ORF',
        back_populates='protein'
    )

    def __init__(self):
        self.residues = []

    @reconstructor
    def init_on_load(self):
        self.residues = [Residue(self, aa, i) for i, aa in enumerate(self.sequence, start=1)]
        if self.orf and self.orf.nucleotides:
            self._link_aa_to_orf_nt()
    
    def __repr__(self):
        return f'{self.orf.transcript}:protein'
    
    @hybrid_property
    def gene(self):
        return self.orf.transcript.gene
    
    @hybrid_property
    def transcript(self):
        return self.orf.transcript

    def _link_aa_to_orf_nt(self):
        aa_sequence = Seq(self.sequence)

        nt_sequence = Seq(self.orf.sequence)
        translation = nt_sequence.translate(to_stop=True)
        aa_match_index = aa_sequence.find(translation)
        if aa_match_index == -1:
            warn(
                f'Could not match amino acid sequence to nucleotide sequence of {self.orf}'
            )
            return
        
        nt_match_index = aa_match_index*3
        nt_list = self.orf.nucleotides[nt_match_index:]
        for i, aa in enumerate(self.residues):
            aa.codon = tuple(nt_list[3*i:3*i + 3])
            for nt in aa.codon:
                nt.residue = aa


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
