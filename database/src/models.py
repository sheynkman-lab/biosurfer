
import enum
from itertools import islice
from warnings import warn

from Bio.Seq import Seq
from inscripta.biocantor.gene.cds import CDSFrame, CDSInterval
from inscripta.biocantor.location.location_impl import (CompoundInterval,
                                                        SingleInterval, Strand)
from inscripta.biocantor.parent import Parent
from inscripta.biocantor.sequence import Alphabet, Sequence
from sqlalchemy import CHAR, Column, Enum, ForeignKey, Integer, String
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.orm import reconstructor, relationship

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
    gene_id = Column(Integer, ForeignKey('gene.id'))
    gene = relationship('Gene', back_populates='transcripts')
    # experiment
    # sample
    exons = relationship(
        'Exon',
        order_by='Exon.start',
        back_populates='transcript',
        uselist=True)
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

    @reconstructor
    def init_on_load(self):
        self.nucleotides = [Nucleotide(self, None, i, nt) for i, nt in enumerate(self.sequence)]
    
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
    def protein(self):
        """Get the "primary" protein produced by this transcript, if it exists."""
        # TODO: implement this
        raise NotImplementedError


class GencodeTranscript(Transcript):
    def __init__(self, *, accession, name, strand, appris, start_nf, end_nf):
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
    # genomic coordinates
    start = Column(Integer)
    stop = Column(Integer)
    # transcript coordinates
    transcript_start = Column(Integer)
    transcript_stop = Column(Integer)
    sequence = Column(String, default='')
    transcript_id = Column(Integer, ForeignKey('transcript.id'))
    transcript = relationship(
        'Transcript', 
        back_populates='exons'
    )
    
    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': 'exon'
    }

    @reconstructor
    def init_on_load(self):
        # TODO: deprecate this after implementing nucleotides property
        self.nucleotides = []
        for i in range(len(self.sequence)):
            nuc_str = self.sequence[i]
            if self.strand is Strand.MINUS:
                coord = self.stop - i
            else:
                coord = self.start + i
            nucleotide = Nucleotide(self, coord, i, nuc_str)
            self.nucleotides.append(nucleotide)

    def __repr__(self) -> str:
        # TODO: change to exon number
        return f'{self.transcript}|{self.start}-{self.stop}'

    @hybrid_property
    def _location(self):
        return SingleInterval(self.transcript_start-1, self.transcript_stop, Strand.PLUS)
    
    # TODO: should pull slice from self.transcript.nucleotides, but Exon.transcript_start and Exon.transcript_stop need to be implemented first
    # @hybrid_property
    # def nucleotides(self):
    #     pass

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
    
    # @hybrid_property
    # def sequence(self):
    #     seq = [nuc.nucleotide for nuc in self.nucleotides]
    #     seq = ''.join(seq)
    #     return seq


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
    def __init__(self, parent, coordinate, position, nucleotide) -> None:
        self.parent = parent
        self.coordinate = coordinate  # genomic coordinate
        self.position = position  # position within parent
        self.nucleotide = nucleotide  # TODO: make this an Enum?
        self.amino_acid = None  # associated AminoAcid, if any
    
    def __repr__(self) -> str:
        return self.nucleotide
    
    @property
    def gene(self):
        if isinstance(self.parent, Gene):
            return self.parent
        return self.parent.gene
    

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
        back_populates='orf'
    )

    def __repr__(self) -> str:
        return f'{self.transcript}|orf:{self.transcript_start}-{self.transcript_stop}'

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
        self.amino_acids = []

    @reconstructor
    def init_on_load(self):
        self.amino_acids = []
        for i, residue in enumerate(self.sequence):
            residue = self.sequence[i]
            amino_acid = AminoAcid(self, residue, i+1)
            self.amino_acids.append(amino_acid)
        if self.orf and self.orf.nucleotides:
            self._link_aa_to_orf_nt()
    
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
        for i, aa in enumerate(self.amino_acids):
            aa.codon = tuple(nt_list[3*i:3*i + 3])
            for nt in aa.codon:
                nt.amino_acid = aa
            # tl = CODON_TABLE[''.join(nt.nucleotide for nt in aa.codon)]
            # assert tl == aa.amino_acid, f'{aa.codon} does not translate to {aa}'


class AminoAcid:
    def __init__(self, protein, amino_acid, position) -> None:
        self.amino_acid = amino_acid  # TODO: make this an Enum?
        self.protein = protein
        self.position = position  # position within protein peptide sequence
        self.codon = (None, None, None)  # 3-tuple of associated Nucleotides; filled in later
    
    def __repr__(self) -> str:
        return f'{self.amino_acid}{self.position}'
