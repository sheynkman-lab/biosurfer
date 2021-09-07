from sqlalchemy import (Column, Integer, ForeignKey, String, Table)
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.orm import (reconstructor, relation, relationship)

from Bio.Seq import Seq

from biosurfer.core.models import (Transcript, Exon, ORF, Protein,
        Nucleotide, Base, AccessionMixin)
from biosurfer.core.constants import Nucleobase, Strand







variant_transcript_association_table = Table('variant_transcript_association', Base.metadata,
    Column('variant_id', ForeignKey('variant.id'), primary_key=True),
    Column('transcript_id', ForeignKey('transcript.accession'), primary_key=True)
)

class Variant(Base):
    """Nucloetide Variant class
    A variant is a modification of the anchor nucleotide
     sequence to the variant nucloetide sequence
    """
    id = Column(Integer, primary_key=True)
    chromosome_id = Column(Integer, ForeignKey('chromosome.name'))
    chromosome_position = Column(Integer)
    reference_nucleotide_sequence = Column(String)
    variant_nucleotide_sequence = Column(String)

    chromosome = relationship('Chromosome')

    def __repr__(self) -> str:
        return f'{self.chromosome}.{self.chromosome_position}'
    variant_transcripts = relationship(
        'VariantTranscript',
        secondary=variant_transcript_association_table,
        back_populates='variants' )


class VariantTranscript(Transcript):
    __mapper_args__ = {
        'polymorphic_identity': 'varianttranscript'
    }
    variants = relationship(
        'Variant',
        secondary=variant_transcript_association_table,
        back_populates='variant_transcripts'

    )


class VariantBuilder:
    def __init__(self, reference_transcript, variants):
        self.variants = variants
        self.reference_transcript = reference_transcript
        self.variant_id = '.'.join([str(variant) for variant in variants])


        self.variant_transcript = self.build_variant_transcript()
        self.variant_exons = self.build_variant_exons()
    
    def build_variant_transcript(self):
        nucleotides = self.reference_transcript.nucleotides.copy()
        variant_nucleotides = []
        for variant in self.variants:
            base_nucleotide = self.reference_transcript.get_nucleotide_from_coordinate(variant.chromosome_position)
            if base_nucleotide is not None:
                variant_nucleotide = Nucleotide(
                    self.reference_transcript, 
                    variant.chromosome_position, 
                    base_nucleotide.position, 
                    Nucleobase(variant.variant_nucleotide_sequence))
                nucleotides[base_nucleotide.position - 1] = variant_nucleotide
                variant_nucleotides.append(variant_nucleotide)
        sequence = ''.join([str(nuc.base) for nuc in nucleotides])
        # self.variant_id = '.'.join([str(nuc) for nuc in variant_nucleotides])
        
        accession = f'{self.reference_transcript.accession}.{self.variant_id}'
        variant_transcript = VariantTranscript(
            accession = accession,
            name=self.reference_transcript.name,
            strand = self.reference_transcript.strand,
            gene = self.reference_transcript.gene,
            sequence = sequence)

        return variant_transcript
        
    def build_variant_exons(self):
        variant_exons = []
        for exon in self.reference_transcript.exons:
            variant_exon = Exon(
                position = exon.position,
                start = exon.start,
                stop = exon.stop,
                transcript_start = exon.transcript_start,
                transcript_stop = exon.transcript_stop,
                transcript = self.variant_transcript,
                accession = f'{exon.accession}.{self.variant_id}'
            )
            variant_exons.append(variant_exon)
        self.variant_transcript.exons = variant_exons

    def build_variant_orfs_proteins(self):

        for orf in self.reference_transcript.orfs:

            if self.reference_transcript.strand is Strand.PLUS:
                nucleotide_sequence = self.variant_transcript.sequence[orf.start-1:]
                nucleotide_sequence = Seq(nucleotide_sequence)
                protein_sequence = nucleotide_sequence.translate(to_stop=True)
            elif self.reference_transcript.strand is Strand.MINUS:
                nucleotide_sequence = self.variant_transcript.sequence[:orf.transcript_stop]
                nucleotide_sequence = Seq(nucleotide_sequence)
                nucleotide_sequence = nucleotide_sequence.reverse_complement()
                protein_sequence = nucleotide_sequence.translate(to_stop=True)
            # orf_stop = len(protien_sequence) + orf_start


    def build_variant_orfs(self):
        pass

# class Variant:
#     """Nucloetide Variant class
#     A variant is a modification of the anchor nucleotide
#      sequence to the variant nucloetide sequence
#     """
#     def __init__(self, chromosome, chromosome_position,
#             reference_nucleotides, variant_nucleotides):
#         self.chromosome = chromosome
#         self.chromosome_position = chromosome_position
#         self.reference_nucleotides = reference_nucleotides
#         self.variant_nucleotides = variant_nucleotides

# class VariantTranscript:
#     def __init__(self, reference_transcript, variants) -> None:
#         self.reference_transcript = reference_transcript
#         self.variants = variants