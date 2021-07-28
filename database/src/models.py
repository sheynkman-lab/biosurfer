
from sqlalchemy import create_engine
from sqlalchemy import Table
from sqlalchemy import Column, Integer, String, Text, ForeignKey, Enum, CHAR
from sqlalchemy.ext import hybrid
from sqlalchemy.orm import declarative_base, relation
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker
from sqlalchemy.orm import reconstructor
import enum
import numpy as np
from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from database import Base
import itertools


# class Strand(enum.Enum):
#     plus = '+'
#     minus = '-'

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
    strand = Column(String)
    transcripts = relationship('Transcript', back_populates='gene', order_by='Transcript.name')
    chromosome = relationship('Chromosome', back_populates='genes')
    def __repr__(self) -> str:
        return self.name


class Transcript(Base):
    __tablename__ = 'transcript'
    id = Column(Integer, primary_key=True)
    accession = Column(String)
    name = Column(String)
    gene_id = Column(Integer, ForeignKey('gene.id'))
    gene = relationship('Gene', back_populates='transcripts')
    # experiment
    # sample
    exons = relationship(
        'Exon',
        order_by='Exon.start',
        back_populates='transcript')
    orfs = relationship(
        'ORF',
        order_by='ORF.start',
        back_populates='transcript',
        uselist=True)

    def __repr__(self) -> str:
        return self.name

    @hybrid_property
    def start(self):
        start = np.inf
        for exon in self.exons:
            if exon.start < start:
                start = exon.start
        return start
    
    @hybrid_property
    def length(self):
        length = sum([exon.length for exon in self.exons])
        return length
    
    @hybrid_property
    def sequence(self):
        sequence = [exon.sequence for exon in self.exons]
        sequence = ''.join(sequence)
        return sequence

    
    @hybrid_property
    def chromosome(self):
        return self.gene.chromosome
    
    @hybrid_property
    def nucleotides(self):
        nucleotides = []
        for exon in self.exons:
            nucleotides = nucleotides + exon.nucleotides
        return nucleotides
    
    @hybrid_property
    def protein(self):
        """Get the "primary" protein produced by this transcript, if it exists."""
        pass
    

class Exon(Base):
    __tablename__ = 'exon'
    id = Column(Integer, primary_key=True)
    accession = Column(String)
    # genomic coordinates
    start = Column(Integer)
    stop = Column(Integer)
    # transcript coordinates
    start_tx = Column(Integer)
    stop_tx = Column(Integer)
    sequence = Column(String, default='')

    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    transcript = relationship(
        'Transcript', 
        back_populates='exons')

    def __init__(self):
        self.nucleotides = []

    @reconstructor
    def init_on_load(self):
        self.nucleotides = []
        for i in range(len(self.sequence)):
            nuc_str = self.sequence[i]
            if self.strand == '-':
                coord = self.stop - i
            else:
                coord = self.start + i
            nucleotide = Nucleotide(self, coord, i, nuc_str)
            self.nucleotides.append(nucleotide)
    
    # nucleotides = relationship(
    #     'Nucleotide',
    #     order_by='Nucleotide.position',
    #     secondary=exon_nucleotide_table,
    #     back_populates='exons')
    def __repr__(self) -> str:
        # TODO: change to exon number
        return f'{self.transcript}|{self.start}-{self.stop}'
    
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
        return self.gene.strand
    
    # @hybrid_property
    # def sequence(self):
    #     seq = [nuc.nucleotide for nuc in self.nucleotides]
    #     seq = ''.join(seq)
    #     return seq


class Nucleotide:
    def __init__(self, exon, coordinate, position, nucleotide) -> None:
        self.exon = exon
        self.coordinate = coordinate  # genomic coordinate
        self.position = position  # position within exon
        self.nucleotide = nucleotide  # TODO: make this an Enum?
    
    def __repr__(self) -> str:
        return self.nucleotide
    
    @property
    def gene(self):
        return self.exon.gene
    

class ORF(Base):
    __tablename__ = 'orf'
    id = Column(Integer, primary_key=True)
    # genomic coordinates
    start = Column(Integer)  
    stop = Column(Integer)  
    # transcript coordinates
    start_tx = Column(Integer)
    stop_tx = Column(Integer)

    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    transcript = relationship(
        'Transcript', 
        back_populates='orfs'
    )
    protein = relationship(
        'Protein',
        back_populates='orf'
    )

    @hybrid_property
    def sequence(self):
        return self.transcript.sequence[self.start_tx - 1:self.stop_tx]
    
    # TODO: uncomment after implementing ORF.start and ORF.stop
    # @reconstructor
    # def init_on_load(self):
    #     self.nucleotides = []
    #     for i in range(len(self.sequence)):
    #         nuc_str = self.sequence[i]
    #         if self.strand == '-':
    #             coord = self.stop - i
    #         else:
    #             coord = self.start + i
    #         nucleotide = Nucleotide(self, coord, i, nuc_str)
    #         self.nucleotides.append(nucleotide)


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
        self.amino_acids = []

    @reconstructor
    def init_on_load(self):
        self.amino_acids = []
        for i in range(len(self.sequence)):
            residue = self.sequence[i]
            nucleotides = []
            amino_acid = AminoAcid(self, residue, i+1)
            self.amino_acids.append(amino_acid)


class AminoAcid():
    def __init__(self, protein, amino_acid, position) -> None:
        self.amino_acid = amino_acid  # TODO: make this an Enum?
        self.protein = protein
        self.position = position  # position within protein_isoform peptide sequence
        self.codon = (None, None, None)
    
    def __repr__(self) -> str:
        return f'{self.amino_acid}{self.position}'
