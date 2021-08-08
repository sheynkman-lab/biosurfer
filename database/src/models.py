from collections import Counter
from operator import attrgetter
from typing import List
from warnings import warn

from Bio.Seq import Seq
from inscripta.biocantor.location.location_impl import (CompoundInterval,
                                                        SingleInterval, Strand)
from sqlalchemy import CHAR, Column, Enum, ForeignKey, Integer, String
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import reconstructor, relationship

from constants import Nucleobase, AminoAcid
from database import Base


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

    def __init__(self):
        self._nucleotides = None

    @reconstructor
    def init_on_load(self):
        self._nucleotides = None

    @hybrid_property
    def nucleotides(self):
        if not self.sequence:
            raise AttributeError(f'{self} has no sequence')
        elif self._nucleotides is None:
            assert sum(exon.stop - exon.start + 1 for exon in self.exons) == len(self.sequence)
            self._nucleotides = []
            i = 0
            for exon in self.exons:
                coords = range(exon.start, exon.stop+1)
                if self.strand is Strand.MINUS:
                    coords = reversed(coords)
                for coord in coords:
                    self._nucleotides.append(Nucleotide(self, coord, i+1, self.sequence[i]))
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
    
    @hybrid_method
    def get_exon_containing_position(self, position: int) -> 'Exon':
        """Given a position (1-based) within the transcript's nucleotide sequence, return the exon containing that position."""
        for exon in self.exons:
            if position in range(exon.transcript_start, exon.transcript_stop + 1):
                return exon
        raise ValueError(f'Position {position} not found in {self}')


class GencodeTranscript(Transcript):
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
    def _location(self):
        return SingleInterval(self.transcript_start-1, self.transcript_stop, Strand.PLUS)


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


class Residue:
    def __init__(self, protein: 'Protein', amino_acid: str, position: int) -> None:
        self.amino_acid = AminoAcid(amino_acid)
        self.protein = protein
        self.position = position  # position within protein peptide sequence
        self.codon = (None, None, None)  # 3-tuple of associated Nucleotides; filled in later
    
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
