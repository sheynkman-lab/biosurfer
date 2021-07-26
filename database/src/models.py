
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
    name = Column(String)
    chormosome_id = Column(Integer,ForeignKey('chromosome.id'))
    strand = Column(String)
    transcripts = relationship('Transcript', back_populates='gene', order_by='Transcript.name')
    chromosome = relationship('Chromosome', back_populates='genes')
    def __repr__(self) -> str:
        return self.name

transcript_exon_association_table = Table('transcript_exon', Base.metadata,
    Column('transcript_id', Integer, ForeignKey('transcript.id')),
    Column('exon_id', Integer, ForeignKey('exon.id'))
)

class Transcript(Base):
    __tablename__ = 'transcript'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    gene_id = Column(Integer, ForeignKey('gene.id'))
    gene = relationship('Gene', back_populates='transcripts')
    exons = relationship(
        'Exon',
        order_by='Exon.start', 
        secondary=transcript_exon_association_table,
        back_populates='transcripts')
    # orf = relationship('ORF', back_populates='transcript', uselist=False)

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
    


# exon_nucleotide_table = Table('exon_nucleotide',  Base.metadata,
#     Column('exon_id', Integer, ForeignKey('exon.id')),
#     Column('nucleotide_id', Integer, ForeignKey('nucleotide.id')))

class Exon(Base):
    __tablename__ = 'exon'
    id = Column(Integer, primary_key=True)
    start = Column(Integer)
    stop = Column(Integer)
    sequence = Column(String, default='')
    transcripts = relationship(
        'Transcript', 
        secondary=transcript_exon_association_table, 
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
        return f'{self.start}-{self.stop}'
    
    @hybrid_property
    def length(self):
        return self.stop - self.start + 1
    
    @hybrid_property
    def gene(self):
        if len(self.transcripts)>0:
            return self.transcripts[0].gene
        return None
    
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
        self.coordinate = coordinate
        self.position = position
        self.nucleotide = nucleotide
    
    def __repr__(self) -> str:
        return self.nucleotide
    
    @property
    def gene(self):
        return self.exon.gene
    

# class ORF(Base):
#     __tablename__ = 'orf'
#     id = Column(Integer, primary_key=True)
#     start = Column(Integer)
#     stop = Column(Integer)
#     transcript_id = Column(Integer,ForeignKey('transcript.id'))
#     protein_isoform_id = Column(Integer,ForeignKey('protein_isoform.id'))
#     transcript = relationship(
#         'Transcript', 
#         back_populates='orfs')

# class Protein(Base):
#     __tablename__= 'protein'
#     id = Column(Integer, primary_key=True)
#     name = Column(String)
#     protein_isoforms  = relationship(
#         'ProteinIsoform', 
#         back_populates='protein')

# class ProteinIsoform(Base):
#     __tablename__ = 'protein_isoform'
#     id = Column(Integer, primary_key=True)
#     protein_isoform_id = Column(Integer,ForeignKey('orf.id'))
#     sequence = Column(String)
#     def __init__(self):
#             self.amino_acids = []

#     @reconstructor
#     def init_on_load(self):
#         self.amino_acids = []
#         for i in range(len(self.sequence)):
#             residue = self.sequence[i]
#             nucleotides = []
#             amino_acid = AminoAcid(self, residue, i)
#             self.amino_acids.append(amino_acid)


# class AminoAcid():
#     def __init__(self, protein_isoform, amino_acid, position) -> None:
#         self.amino_acid = amino_acid
#         self.protein_isoform = protein_isoform
#         self.position = position


