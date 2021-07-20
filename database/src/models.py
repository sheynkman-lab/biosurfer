
from sqlalchemy import create_engine
from sqlalchemy import Table
from sqlalchemy import Column, Integer, String, Text, ForeignKey, Enum, CHAR
from sqlalchemy.ext import hybrid
from sqlalchemy.orm import declarative_base, relation
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker
import enum
import numpy as np
from sqlalchemy.ext.hybrid import hybrid_property, hybrid_method
from database import Base


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
    orf = relationship('ORF', back_populates='transcript', uselist=False)

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
    


exon_nucleotide_table = Table('exon_nucleotide',  Base.metadata,
    Column('exon_id', Integer, ForeignKey('exon.id')),
    Column('nucleotide_id', Integer, ForeignKey('nucleotide.id')))

class Exon(Base):
    __tablename__ = 'exon'
    id = Column(Integer, primary_key=True)
    start = Column(Integer)
    stop = Column(Integer)
#   sequence...
    transcripts = relationship(
        'Transcript', 
        secondary=transcript_exon_association_table, 
        back_populates='exons')
    
    nucleotides = relationship(
        'Nucleotide',
        order_by='Nucleotide.coordinate',
        secondary=exon_nucleotide_table,
        back_populates='exons')
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
    def sequence(self):
        seq = [nuc.nucleotide for nuc in self.nucleotides]
        seq = ''.join(seq)
        return seq

orf_cds_association_table = Table('orf_cds', Base.metadata,
    Column('orf_id', Integer, ForeignKey('orf.id')),
    Column('cds_id', Integer, ForeignKey('cds.id'))
)

class Nucleotide(Base):
    __tablename__ = 'nucleotide'
    id = Column(Integer, primary_key=True)
    coordinate = Column(Integer)
    nucleotide = Column(CHAR)
    exons = relationship(
        'Exon', 
        secondary=exon_nucleotide_table, 
        back_populates='nucleotides')

    def __repr__(self) -> str:
        return self.nucleotide

class ORF(Base):
    __tablename__ = 'orf'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    transcript_id = Column(Integer, ForeignKey('transcript.id'))

    transcript = relationship('Transcript', back_populates='orf')
    cds = relationship(
        'CDS',
        order_by='CDS.start', 
        secondary=orf_cds_association_table,
        back_populates='orfs')
    def __repr__(self) -> str:
        return self.name


class CDS(Base):
    __tablename__ = 'cds'
    id = Column(Integer, primary_key=True)
    start = Column(Integer)
    stop = Column(Integer)

    orfs = relationship(
        'ORF', 
        secondary=orf_cds_association_table, 
        back_populates='cds')
    def __repr__(self) -> str:
        return f'{self.start}-{self.stop}'
