from typing import Dict, Iterable, List, Optional, Tuple
from functools import cached_property

from sqlalchemy import (Boolean, Column, Enum, ForeignKey, Integer, String, Table,
                        create_engine)
from sqlalchemy.orm import (reconstructor, relationship, scoped_session,
                            sessionmaker)
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property

from biosurfer.core.database import Base
from biosurfer.core.constants import APPRIS, SQANTI, AminoAcid, Nucleobase, Strand
from biosurfer.core.models import Transcript, Exon, Protein, Nucleotide, Residue


class VariantTranscriptLink(Base):
    variant_id = Column(Integer, ForeignKey('variant.id'), primary_key=True)
    transcript_id = Column(Integer, ForeignKey('transcript.accession'), primary_key=True)


class Variant(Base):
    """Nucloetide Variant class
    A variant is a modification of the anchor nucleotide
     sequence to the variant nucloetide sequence
    """
    id = Column(Integer, primary_key=True)
    chromosome_id = Column(Integer, ForeignKey('chromosome.name'))
    position = Column(Integer)
    reference_sequence = Column(String)
    variant_sequence = Column(String)
    # quality_score = Column(Float)


    chromosome = relationship('Chromosome')

    def __repr__(self) -> str:
        return f'{self.chromosome}.{self.position}[{self.reference_sequence}->{self.variant_sequence}]'
    variant_transcripts = relationship(
        'VariantTranscript',
        secondary='varianttranscriptlink',
        back_populates='variants')
    
class VariantTranscript(Transcript):
    reference_transcript_id = Column(Integer, ForeignKey('transcript.accession'))
    variants = relationship(
        'Variant',
        secondary='varianttranscriptlink',
        back_populates='variant_transcripts'
    )
    __mapper_args__ = {
        'polymorphic_identity': 'varianttranscript'
    }
    @hybrid_property
    def reference_transcript(self):
        return Transcript.from_accession(self.reference_transcript_id)
    
    @hybrid_property
    def variant_positions(self):
        return [var.position for var in self.variants]

    
    @cached_property
    def nucleotides(self):
        if not self.sequence:
            raise AttributeError(f'{self.name} has no sequence')
        nucleotides = []
        self._nucleotide_mapping: Dict[int, 'Nucleotide'] = dict()
        i = 0
        for exon in self.exons:
            coords = range(exon.start, exon.stop+1)
            if self.strand is Strand.MINUS:
                coords = reversed(coords)
            for coord in coords:
                if coord in self.variant_positions:
                    



                nt = Nucleotide(self, coord, i+1, self.sequence[i])
                nucleotides.append(nt)
                self._nucleotide_mapping[coord] = nt
                i += 1
        return nucleotides

    def get_variant_nucleotides(self):
        pass
    def get_variant_nucleotide(self, variant):
        if variant in self.variants:
            pass
        else :
            raise Exception('variant not in transcript')
    
class VariantExon(Exon):
    __mapper_args__ = {
        'polymorphic_identity': 'variantexon'
    }
    def get_variants(self):
        pass
    def get_variant_nucleotides(self):
        pass
    def get_variant_nucleotide(self, variant):
        pass

class VariantProtein(Protein):
    __mapper_args__ = {
        'polymorphic_identity': 'variantprotein'
    }
    reference_protein_id = Column(Integer, ForeignKey('protein.accession'))

    def get_variant_residues(self):
        pass
    def get_variant_peptides(self):
        pass
    def get_synonymous_residues(self):
        pass
    def get_non_synonymous_residues(self):
        pass

class VariantNucleotide(Nucleotide):
    def __init__(self, parent, coordinate: int, position: int, base: str,  variant: Variant) -> None:
        super().__init__(parent, coordinate, position, base )
        self.variant = variant

class VariantResidue(Residue):
    def __init__(self, protein: 'Protein', amino_acid: str, position: int, variant: Variant) -> None:
        super().__init__(protein, amino_acid, position)
        self.variant = Variant