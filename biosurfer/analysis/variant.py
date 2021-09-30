#%%
import os
from tqdm import tqdm
from sqlalchemy.sql.expression import desc
from biosurfer.core.models import Protein, Transcript, VariantTranscript, VariantTranscriptLink
from biosurfer.core.variant import (
    get_transcript_accessions_with_variants,
    get_possible_variants,
    get_possible_exonic_variants,
    VariantBuilder
)
from biosurfer.core.database import db_session
from biosurfer.core.populate_database import load_gencode_gtf, load_transcript_fasta, load_translation_fasta, load_variants, FastaHeaderFields
import logging

results_dir = '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/biosurfer/biosurfer/analysis/results'

logging.basicConfig(filename=os.path.join(results_dir,'variant_logging.log'), encoding='utf-8', level=logging.DEBUG)
#%%

def get_ids_from_gencode_fasta(header: str):
    fields = [field for field in header.split('|') if field and not field.startswith(('UTR', 'CDS'))]
    transcript_id = next((field for field in fields if field.startswith('ENST')))
    protein_id = next((field for field in fields if field.startswith('ENSP')), None)
    return FastaHeaderFields(transcript_id, protein_id)

def skip_par_y(header: str):  # these have duplicate ENSEMBL accessions and that makes SQLAlchemy very sad
    return 'PAR_Y' in header

def populate_database():
    # Populate Database
    gtf_file = '../../data/gencode.v38.basic.annotation.gtf'
    transcript_file = '../../data/gencode.v38.pc_transcripts.fa'
    translation_file = '../../data/gencode.v38.pc_translations.fa'
    variants_file = '../../data/both_isoseq_shortreads.vcf'
    load_gencode_gtf(gtf_file, overwrite=True)
    load_transcript_fasta(transcript_file, get_ids_from_gencode_fasta, skip_par_y)
    load_translation_fasta(translation_file, get_ids_from_gencode_fasta, skip_par_y)

    load_variants(variants_file)

def build_transcript_variants():
    transcripts_with_possible_variants = get_transcript_accessions_with_variants()
    print(len(transcripts_with_possible_variants))
    i = 0
    t = tqdm(
        transcripts_with_possible_variants,
        desc = 'Finding unique variant peptides...',
        total = len(transcripts_with_possible_variants),
        unit = 'transcripts'
    )
    for i, accession in enumerate(t):
        accession = accession[0]
        try:
            transcript = Transcript.from_accession(accession)
            possible_variants = get_possible_exonic_variants(transcript)
            if len(possible_variants)>0:
                variant_builder = VariantBuilder(transcript,possible_variants)
            # else:
            #     print(f'No Exonic Variants\t{transcript}')
        except:
            logging.error(f'Bad Transcript\t{transcript.accession}')


def find_variant_peptides(variant_transcript):
    """[summary]

    Args:
        variant_transcript ([type]): [description]

    Returns:
        [type]: [description]
    """
    digest = ProteinDigest(missed_cleavages=0)
    peptides = digest.digest_protein(variant_transcript.protein)
    variant_peptides = []
    for variant in variant_transcript.variants:
        variant_nucleotide = variant_transcript.get_nucleotide_from_coordinate(variant.position)
        variant_residue = variant_nucleotide.residue
        for peptide in peptides:
            if variant_residue in peptide.residues:
                variant_peptides.append({
                    'transcriipt': variant_transcript,
                    'variant': variant,
                    'peptide': peptide
                })
    return variant_peptides


def filter_unique_variant_peptides(variant_transcript, variant_peptides):
    unique_peptides = []
    digest = ProteinDigest(missed_cleavages=0)
    ref_peptides = digest.digest_protein(variant_transcript.reference_transcript.protein)
    ref_pep_str = [str(ref_pep) for ref_pep in ref_peptides]
    for var_pep in variant_peptides:
        var_pep_str = str(var_pep['peptide'])
        if var_pep_str not in ref_pep_str:
            unique_peptides.append(var_pep)
    return unique_peptides


def find_all_unique_variant_peptides():
    variant_transcripts = db_session.query(VariantTranscript).all()
    unique_peptides = []
    variant_peptides = []
    t = tqdm(
        variant_transcripts,
        desc = 'Finding unique variant peptides...',
        total = len(variant_transcripts),
        unit = 'transcripts'
    )
    for i, vt in enumerate(t):
        try:
            variant_peps = find_variant_peptides(vt)
            unique_peps = filter_unique_variant_peptides(vt, variant_peps)
            variant_peptides.append(variant_peps)
            unique_peptides.append(unique_peps)
        except:
            logging.error(f'Error finding variant peptides. Variant Transcript: \t{vt.accession}')
    return variant_peptides, unique_peptides

def build_variant_peptide_database(unique_pepties, filename):
        t = tqdm(
        unique_peptides,
        desc = 'Writing unique variant peptides...',
        total = len(unique_peptides),
        unit = 'peptides'
        )
        with open(filename, 'w') as ofile:
            for i, unique_pep in enumerate(t):
                if len(unique_pep.residues)>7:
                    ofile.write(
                        f'>vt|{unique_pep["transcript"].reference_transcriipt.accession}.{unique_pep["variant"].variant_id}\n' + \
                            f'{str(unique_pep["peptide"])}\n'
                    )


def build_variant_protein_databse(unique_peptides, filename):
    unique_transcripts = []
    for unique_pep in unique_peptides:
        unique_transcripts.append(unique_pep['transcript'])
    unique_transcripts = set(unique_transcripts)
    with open(filename) as ofile:
        for ut in unique_transcripts:
            ofile.write(f'>vt|{ut.reference_transcript.accession}\n{ut.protein.sequence}\n')

#%%
# populate_database()
# build_transcript_variants() 
#%%           
variant_peptides, unique_peptides = find_all_unique_variant_peptides()

build_variant_peptide_database(unique_peptides, os.path.join(results_dir, 'variant_peptides.fasta'))
build_variant_protein_databse(unique_peptides, os.path.join(results_dir, 'variant_proteins.fasta'))



# %%



def build_peptide_metadata_table(variant_peptides, category, filename):
    header = [
        'Accession', 'Category', 'Protein Change', 
        'Reference Peptide', 'Variant Peptide', 
        'Protein Accession', 'Protein Position', 'Strand',
        'Chromosome', 'Coordinate', 
        'Reference Nucleotide', 'Variant Nucloetide', 
        'Reference Codon', 'Variant Codon',
    ]

from biosurfer.core.helpers import StringEnum
class VariantCategory(StringEnum):
    INTRON='intronic'
    FIVE_PRIME_UTR='five-prime-utr'
    THREE_PRIME_UTR='three-prime-utr'
    SYNONYMOUS_CODING='synonymous coding'
    NON_SYNONYMOUS_CODING='non-synonymous coding'



def categorize_variant_in_transcript(variant, transcript):

def build_variant_metadata_table():
    header = [
        'Variant', 'Transcript','Category',
    ]