from functools import cached_property
from operator import attrgetter
from typing import Dict, Iterable, List, Optional, Tuple, Type
from warnings import warn

from Bio.Seq import Seq
from biosurfer.core.constants import APPRIS, SQANTI, Strand
from biosurfer.core.helpers import BisectDict
from biosurfer.core.models.base import AccessionMixin, Base, NameMixin, TablenameMixin
from biosurfer.core.models.nonpersistent import *
from more_itertools import only
from sqlalchemy import Boolean, Column, Enum, ForeignKey, Integer, String, func
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.ext.orderinglist import ordering_list
from sqlalchemy.orm import relationship


# TODO: replace with Enum?
class Chromosome(Base, TablenameMixin, NameMixin):
    name = Column(String, primary_key=True)
    genes = relationship('Gene', back_populates='chromosome')

    def __repr__(self) -> str:
        return self.name


class Gene(Base, TablenameMixin, NameMixin, AccessionMixin):
    strand = Column(Enum(Strand))
    chromosome_id = Column(String, ForeignKey('chromosome.name'))
    chromosome = relationship('Chromosome', back_populates='genes')
    transcripts = relationship(
        'Transcript',
        back_populates = 'gene',
        order_by = 'Transcript.name',
        lazy = 'selectin'  # always load transcripts along with gene
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
    


class Transcript(Base, TablenameMixin, NameMixin, AccessionMixin):
    strand = Column(Enum(Strand))
    type = Column(String)
    sequence = Column(String)
    gene_id = Column(String, ForeignKey('gene.accession'))
    gene = relationship('Gene', back_populates='transcripts')
    exons = relationship(
        'Exon',
        order_by = 'Exon.transcript_start',
        collection_class = ordering_list('position', count_from=1),
        back_populates = 'transcript',
        uselist = True,
        lazy = 'selectin'  # always load exons along with transcript
    )
    orfs = relationship(
        'ORF',
        order_by = 'ORF.transcript_start',
        back_populates = 'transcript',
        uselist = True
    )

    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': 'transcript'
    }

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
            chr = self.gene.chromosome_id
            if self.strand is Strand.MINUS:
                up_exon_stop = up_exon.start
                down_exon_start = down_exon.stop
            else:
                up_exon_stop = up_exon.stop
                down_exon_start = down_exon.start
            donor = Position(chr, self.strand, up_exon_stop) + 1
            acceptor = Position(chr, self.strand, down_exon_start) - 1
            junction = Junction.from_splice_sites(donor, acceptor)
            mapping[junction] = (up_exon, down_exon)
        return mapping
    
    @cached_property
    def nucleotides(self):
        if not self.sequence:
            raise AttributeError(f'{self.name} has no sequence')
        nucleotides = []
        self._nucleotide_mapping: Dict[int, 'Nucleotide'] = dict()
        if self.exons:
            i = 0
            for exon in self.exons:
                coords = range(exon.start, exon.stop+1)
                if self.strand is Strand.MINUS:
                    coords = reversed(coords)
                for coord in coords:
                    nt = Nucleotide(self, coord, i+1)
                    nucleotides.append(nt)
                    self._nucleotide_mapping[coord] = nt
                    i += 1
        else:
            nucleotides = [Nucleotide(self, None, i+1) for i in range(len(self.sequence))]
        return nucleotides

    def __repr__(self) -> str:
        return self.name

    # FIXME: make this work in SQL queries
    @property
    def start(self) -> int:
        return self.exons[0].stop if self.strand is Strand.MINUS else self.exons[0].start
    
    @property
    def stop(self) -> int:
        return self.exons[-1].start if self.strand is Strand.MINUS else self.exons[-1].stop
    
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
        if not hasattr(self, '_nucleotide_mapping'):
            self.nucleotides
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

    def get_genome_coord_from_transcript_coord(self, tx_coord: int) -> Position:
        try:
            nt = self.nucleotides[tx_coord]
            return Position(self.gene.chromosome_id, self.strand, nt.coordinate)
        except IndexError:
            pos = Position(self.gene.chromosome_id, self.strand, self.nucleotides[-1].coordinate)
            return pos + (tx_coord - self.length + 1)
    
    def get_transcript_coord_from_genome_coord(self, gn_coord: Position) -> int:
        if self.strand is not gn_coord.strand:
            raise ValueError(f'{gn_coord} is different strand from {self}')
        elif self.gene.chromosome_id != gn_coord.chromosome:
            raise ValueError(f'{gn_coord} is different chromosome from {self}')
        nt = self.get_nucleotide_from_coordinate(gn_coord.coordinate)
        if nt:
            return nt.position - 1
        else:
            raise KeyError
        

class GencodeTranscript(Transcript):
    __tablename__ = None
    appris = Column(Enum(APPRIS, values_callable=lambda x: [str(m.value) for m in x]))
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
    __tablename__ = None
    sqanti = Column(Enum(SQANTI, values_callable=lambda x: [str(m.value) for m in x]))
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


class Exon(Base, TablenameMixin, AccessionMixin):
    # type = Column(String)
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
        back_populates = 'exons'
    )

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


class ORF(Base, TablenameMixin):
    # genomic coordinates
    # start = Column(Integer)  
    # stop = Column(Integer)  
    # transcript coordinates
    transcript_start = Column(Integer)
    transcript_stop = Column(Integer)
    has_stop_codon = Column(Boolean)
    position = Column(Integer, primary_key=True)
    transcript_id = Column(String, ForeignKey('transcript.accession'), primary_key=True)
    transcript = relationship(
        'Transcript', 
        back_populates = 'orfs',
        lazy = 'joined'
    )
    protein_id = Column(String, ForeignKey('protein.accession'))
    protein = relationship(
        'Protein',
        back_populates = 'orf',
        uselist = False,
        lazy = 'joined'
    )

    @property
    def _first_exon_index(self):
        return self.transcript.get_exon_index_containing_position(self.transcript_start)

    @property
    def _last_exon_index(self):
        try:
            return self.transcript.get_exon_index_containing_position(self.transcript_stop)
        except KeyError as e:
            warn(
                f'KeyError: {e} when getting ORF._last_exon_index for {self}'
            )
            return len(self.transcript.exons) - 1

    @cached_property
    def utr5(self):
        utr5_boundary_exon_index = self._first_exon_index
        if self.transcript.exons[utr5_boundary_exon_index].transcript_start == self.transcript_start:
            utr5_boundary_exon_index -= 1
        if utr5_boundary_exon_index >= 0:
            return FivePrimeUTR(self, utr5_boundary_exon_index)
        else:
            return None
    
    @cached_property
    def utr3(self):
        utr3_boundary_exon_index = self._last_exon_index
        if self.transcript.exons[utr3_boundary_exon_index].transcript_stop == self.transcript_stop:
            utr3_boundary_exon_index += 1
        if utr3_boundary_exon_index < len(self.transcript.exons):
            return ThreePrimeUTR(self, utr3_boundary_exon_index)
        else:
            return None

    @cached_property
    def nmd(self):
        # ORFs with stop codons at least 50 bp upstream of the last splice site in the mature transcript
        # (i.e. the beginning of the last exon) are considered candidates for nonsense-mediated decay (NMD)
        last_junction = self.transcript.exons[-1].transcript_start
        return last_junction - self.transcript_stop >= 50

    def __repr__(self) -> str:
        return f'{self.transcript}:orf({self.transcript_start}-{self.transcript_stop})'

    # FIXME: make this work in SQL queries
    @property
    def start(self):
        if self.transcript.strand is Strand.PLUS:
            return self.nucleotides[0].coordinate
        elif self.transcript.strand is Strand.MINUS:
            return self.nucleotides[-1].coordinate

    @property
    def stop(self):
        if self.transcript.strand is Strand.PLUS:
            return self.nucleotides[-1].coordinate
        elif self.transcript.strand is Strand.MINUS:
            return self.nucleotides[0].coordinate
    
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


class Protein(Base, TablenameMixin, AccessionMixin):
    sequence = Column(String)
    orf = relationship(
        'ORF',
        back_populates = 'protein',
        uselist = False,
        lazy = 'joined'
    )
    features = relationship(
        'ProteinFeature',
        back_populates = 'protein',
        uselist = True,
        order_by = 'ProteinFeature.protein_start'
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
        # if not self.sequence.endswith('*'):
        #     self.sequence += '*'
        _residues = [Residue(self, i+1) for i in range(len(self.sequence))]
        if self.orf:
            self.orf._link_aa_to_nt(_residues)
        return _residues

    def __repr__(self):
        if self.orf:
            return f'{self.orf.transcript}:protein'
        else:
            return self.accession
    
    @property
    def gene(self):
        return self.orf.transcript.gene
    
    @property
    def transcript(self) -> 'Transcript':
        return self.orf.transcript
    
    @hybrid_property
    def length(self):
        return len(self.sequence)
    
    def get_protein_coord_from_transcript_coord(self, transcript_coord: int):
        return (transcript_coord - self.orf.transcript_start + 1) // 3
    
    def get_transcript_coord_from_protein_coord(self, protein_coord: int):
        return protein_coord*3 + self.orf.transcript_start - 1
