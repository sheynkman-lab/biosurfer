
from sqlalchemy import create_engine
from sqlalchemy import Table
from sqlalchemy import Column, Integer, String, Text, ForeignKey, Enum, CHAR
from sqlalchemy.orm import declarative_base, relation
from sqlalchemy.orm import relationship
from sqlalchemy.orm import sessionmaker
import enum

# create an engine
engine = create_engine('sqlite:///:memory:', echo=False)
# create a configured "Session" class
Session = sessionmaker(bind=engine)
# create a Session
session = Session()
Base = declarative_base()

class Strand(enum.Enum):
    plus = '+'
    minus = '-'

class Chromosome(Base):
    __tablename__ = 'chromosome'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    genes = relationship('Gene', back_populates='chromosome')


class Gene(Base):
    __tablename__ = 'gene'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    chormosome_id = ForeignKey('chromosome.id')
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
    __tablename__ = 'orf'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    gene_id = Column(Integer, ForeignKey('gene.id'))

    gene = relationship('Gene', back_populates='transcripts')
    exons = relationship(
        'Exon',
        order_by='Exon.start', 
        secondary=transcript_exon_association_table,
        back_populates='orf')
    orf = relationship('ORF', back_populates='transcript', uselist=False)
    def __repr__(self) -> str:
        return self.name

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
    def __repr__(self) -> str:
        return f'{self.start}-{self.stop}'

orf_cds_association_table = Table('orf_cds', Base.metadata,
    Column('orf_id', Integer, ForeignKey('orf.id')),
    Column('cds_id', Integer, ForeignKey('cds.id'))
)

class ORF(Base):
    __tablename__ = 'orf'
    id = Column(Integer, primary_key=True)
    name = Column(String)
    transcript_id = Column(Integer, ForeignKey('gene.id'))

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
