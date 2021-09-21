from sqlalchemy import select
from sqlalchemy import (Column, Integer, ForeignKey, String, Table, Float)
from sqlalchemy.ext.hybrid import hybrid_method, hybrid_property
from sqlalchemy.orm import (reconstructor, relation, relationship)

from Bio.Seq import Seq

from biosurfer.core.models import (Transcript, Exon, ORF, Protein, Gene,
        Nucleotide, Base, AccessionMixin, 
        Variant, VariantTranscript, db_session)
from biosurfer.core.constants import Nucleobase, Strand


def get_possible_variant_transcripts(variant:Variant):
    possible_transcripts = (
        db_session.query(Transcript)
        .join(Gene)
        .filter(Gene.chromosome_id == variant.chromosome_id)
        .filter(Transcript.start >= variant.position)
        .filter(Transcript.stop <= variant.position)
        ).all()
    return possible_transcripts

        
def get_possible_variants(transcript:Transcript):
    possible_variants = (
        db_session.query(Variant)
        .filter(Variant.chromosome_id == transcript.chromosome.name)
        .filter(Variant.position.between(transcript.start, transcript.stop))
        .all()
    )
    return possible_variants
        # statement = select(Variant)
        # statement = statement.filter(Variant.chromosome_id == transcript.chromosome.name)
        
        # statement = statement.filter(Variant.position.between(transcript.start, transcript.stop))


def get_transcript_accessions_with_variants():
    possible_transcript_variants = (
        db_session.query(Transcript.accession)
        .join(Gene)
        .filter(Variant.chromosome_id == Gene.chromosome_id)
        .filter(Variant.position.between(Transcript.start, Transcript.stop))
        .group_by(Transcript.accession)
    ).all()
    return possible_transcript_variants


class VariantBuilder:
    def __init__(self, reference_transcript, variants):
        self.variants = variants
        self.reference_transcript = reference_transcript
        self.variant_id = '.'.join([str(variant) for variant in variants])
        if self.reference_transcript.sequence is None:
            return
        self.variant_transcript = self.build_variant_transcript()
        self.variant_exons = self.build_variant_exons()
        self.variant_orfs, self.variant_proteins = self.build_variant_orfs_and_proteins()
        self.update_database()

    
    def update_database(self):
        if VariantTranscript.from_accession(self.variant_transcript['accession']) is not None:
            return
        db_session.bulk_insert_mappings(VariantTranscript, [self.variant_transcript])
        db_session.bulk_insert_mappings(Exon, self.variant_exons)
        db_session.bulk_insert_mappings(ORF, self.variant_orfs)
        db_session.bulk_insert_mappings(Protein, self.variant_proteins)
        db_session.commit()
    
    def extract_variant_data(self):
        self.variant_transcript = VariantTranscript.from_accession(self.variant_transcript['accession'])
        


    def build_variant_transcript(self):
        nucleotides = self.reference_transcript.nucleotides.copy()
        variant_nucleotides = []
        for variant in self.variants:
            base_nucleotide = self.reference_transcript.get_nucleotide_from_coordinate(variant.position)
            if base_nucleotide is not None:
                variant_nucleotide = Nucleotide(
                    self.reference_transcript, 
                    variant.position, 
                    base_nucleotide.position, 
                    Nucleobase(variant.variant_sequence))
                nucleotides[base_nucleotide.position - 1] = variant_nucleotide
                variant_nucleotides.append(variant_nucleotide)
        sequence = ''.join([str(nuc.base) for nuc in nucleotides])
        # self.variant_id = '.'.join([str(nuc) for nuc in variant_nucleotides])
        
        accession = f'{self.reference_transcript.accession}.{self.variant_id}'
        variant_transcript = {
            'accession': accession,
            'name': self.reference_transcript.name,
            'strand' : self.reference_transcript.strand,
            'gene_id' : self.reference_transcript.gene.accession,
            'sequence' :  sequence,
            'reference_transcript_id': self.reference_transcript.accession,
            'type': 'varianttranscript'
        }
        return variant_transcript
        
    def build_variant_exons(self):
        variant_exons = []
        for exon in self.reference_transcript.exons:
            variant_exon = {
                'position' : exon.position,
                'start' : exon.start,
                'stop' : exon.stop,
                'transcript_start' : exon.transcript_start,
                'transcript_stop' : exon.transcript_stop,
                'transcript_id' : self.variant_transcript['accession'],
                'accession' : f'{exon.accession}.{self.variant_id}',
                'type' : 'exon'
            }

            # variant_exon = Exon(
            #     position = exon.position,
            #     start = exon.start,
            #     stop = exon.stop,
            #     transcript_start = exon.transcript_start,
            #     transcript_stop = exon.transcript_stop,
            #     transcript = self.variant_transcript,
            #     accession = f'{exon.accession}.{self.variant_id}'
            # )
            variant_exons.append(variant_exon)
        return variant_exons

    def build_variant_orfs_and_proteins(self):
        variant_orfs = []
        variant_proteins = []

        for orf in self.reference_transcript.orfs:
            transcript_start = orf.transcript_start
            if self.reference_transcript.strand is Strand.PLUS or self.reference_transcript.strand is Strand.MINUS:
                nucleotide_sequence = self.variant_transcript['sequence'][transcript_start-1:]
                nucleotide_sequence = Seq(nucleotide_sequence)
                protein_sequence = nucleotide_sequence.translate(to_stop=True)
                transcript_stop = len(protein_sequence) + transcript_start
            
            # currently appears we do not need to do reverse complement of transcript sequence
            # elif self.reference_transcript.strand is Strand.MINUS:
            #     nucleotide_sequence = self.variant_transcript['sequence'][:orf.transcript_stop]
            #     nucleotide_sequence = Seq(nucleotide_sequence)
            #     nucleotide_sequence = nucleotide_sequence.reverse_complement()
            #     protein_sequence = nucleotide_sequence.translate(to_stop=True)
            #     transcript_stop = orf.transcript_start - len(protein_sequence)
            # orf_stop = len(protein_sequence) + transcript_start
            new_protein = {
                'accession': self.variant_transcript['accession'],
                'sequence': str(protein_sequence)
            }
            new_orf = {
                'transcript_id':self.variant_transcript['accession'],
                'transcript_start':transcript_start,
                'transcript_stop':transcript_stop,
                'protein_id':new_protein['accession']
            }
            variant_orfs.append(new_orf)
            variant_proteins.append(new_protein)
            
        return variant_orfs, variant_proteins
        


            
        




# class Variant:
#     """Nucloetide Variant class
#     A variant is a modification of the anchor nucleotide
#      sequence to the variant nucloetide sequence
#     """
#     def __init__(self, chromosome, position,
#             reference_nucleotides, variant_nucleotides):
#         self.chromosome = chromosome
#         self.position = position
#         self.reference_nucleotides = reference_nucleotides
#         self.variant_nucleotides = variant_nucleotides

# class VariantTranscript:
#     def __init__(self, reference_transcript, variants) -> None:
#         self.reference_transcript = reference_transcript
#         self.variants = variants