#%%
import enum
import os
import logging
import pickle
import unicodedata
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
from biosurfer.core.peptides import ProteinDigest
results_dir = '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/biosurfer/biosurfer/analysis/results'

logging.basicConfig(filename=os.path.join(results_dir,'variant_logging_both_isoseq_short.log'), encoding='utf-8', level=logging.DEBUG)
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
        desc = 'Building variant transcript database...',
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
    reference_transcript = variant_transcript.reference_transcript
    digest = ProteinDigest(missed_cleavages=0)
    var_transcript_peptides = digest.digest_protein(variant_transcript.protein)
    ref_transcript_peptides = digest.digest_protein(reference_transcript.protein)
    variant_peptides = []
    for variant in variant_transcript.variants:
        variant_nucleotide = variant_transcript.get_nucleotide_from_coordinate(variant.position)
        variant_residue = variant_nucleotide.residue

        reference_nucleotide = reference_transcript.get_nucleotide_from_coordinate(variant.position)
        reference_residue = reference_nucleotide.residue

        variant_peptide_dict = {
            'variant_transcript': variant_transcript,
            'reference_transcript': reference_transcript,
            'variant': variant,
            'reference_nucleotide': reference_nucleotide,
            'variant_nucleotide': variant_nucleotide,
        }
        for peptide in var_transcript_peptides:
            if variant_residue in peptide.residues:
                variant_peptide_dict['variant_peptide'] = peptide
                variant_peptide_dict['variant_residue'] = variant_residue
                # variant_peptide_dict['variant_nucleotide'] = variant_nucleotide
                break

        for peptide in ref_transcript_peptides:
            if reference_residue in peptide.residues:
                variant_peptide_dict['reference_peptide'] = peptide
                variant_peptide_dict['reference_residue'] = reference_residue
                # variant_peptide_dict['reference_nucleotide'] = reference_nucleotide
                break
        # if 'variant_peptide' in variant_peptide_dict.keys():
        #     variant_peptides.append(variant_peptide_dict)
        
    return variant_peptides


def filter_unique_variant_peptides(variant_peptides):
    unique_peptides = []
    for var_pep in variant_peptides:
        if 'variant_peptide' in var_pep.keys() and 'reference_peptide' in var_pep.keys():
            if str(var_pep['variant_peptide']) != str(var_pep['reference_peptide']):
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
            variant_peptides = variant_peptides + variant_peps
            unique_peps = filter_unique_variant_peptides(variant_peps)
            unique_peptides = unique_peptides + unique_peps
        except:
            logging.error(f'Error finding variant peptides. Variant Transcript: \t{vt.accession}')
            db_session.rollback()
    return variant_peptides, unique_peptides
#%%
def build_variant_peptide_database(unique_peptides, filename):
        t = tqdm(
        unique_peptides,
        desc = 'Writing unique variant peptides...',
        total = len(unique_peptides),
        unit = 'peptides'
        )
        with open(filename, 'w') as ofile:
            for i, unique_pep in enumerate(t):
                if len(unique_pep["variant_peptide"].residues)>7:
                    ofile.write(
                        f'>vt|{unique_pep["reference_transcript"].accession}.{unique_pep["variant"].id}\n' + \
                            f'{str(unique_pep["variant_peptide"])}\n'
                    )


def build_variant_protein_databse(unique_peptides, filename):
    unique_transcripts = []
    for unique_pep in unique_peptides:
        unique_transcripts.append(unique_pep['variant_transcript'])
    unique_transcripts = set(unique_transcripts)
    with open(filename, 'w') as ofile:
        for ut in unique_transcripts:
            ofile.write(f'>vt|{ut.reference_transcript.accession}\n{ut.protein.sequence}\n')

#%%
populate_database()
build_transcript_variants() 
# #%%           
variant_peptides, unique_peptides = find_all_unique_variant_peptides()

#%%
build_variant_peptide_database(unique_peptides, os.path.join(results_dir, 'variant_peptides_isoseq_only.fasta'))
build_variant_protein_databse(unique_peptides, os.path.join(results_dir, 'variant_proteins_isoseq_only.fasta'))



pickle.dump(
    unique_peptides,
    open(os.path.join(results_dir, "unique_peptides_isoseq_only.p"), "wb")
)

#%%
pickle.dump(
    variant_peptides,
    open(os.path.join(results_dir, "variant_peptides_isoseq_only.p"), "wb")
)
#%%

def get_peptide_string(peptide, small_residue):
    peptide_string = ''
    for res in peptide.residues:
        aa = str(res.amino_acid)
        if res == small_residue:
            peptide_string = peptide_string + aa.lower()
        else:
            peptide_string = peptide_string + aa
    return peptide_string

def get_codon_string(residue, nucleotide):
    codon_string = ''
    for nuc in residue.codon:
        base = str(nuc.base)
        if nuc == nucleotide:
            codon_string = codon_string + base.lower()
        else: 
            codon_string = codon_string + base
    return codon_string

def get_codon_position(residue, nucleotide):
    if nucleotide in residue.codon:
        return residue.codon.index(nucleotide) + 1
    return -1

def build_peptide_metadata_table(variant_peptides, database, filename):
    header = [
        'Accession', 'Category', 'Protein Change', 
        'Reference Peptide', 'Variant Peptide', 
        'Protein Accession', 'Protein Position', 'Strand',
        'Chromosome', 'Coordinate', 
        'Reference Nucleotide', 'Variant Nucloetide', 
        'Reference Codon', 'Variant Codon',
    ]
    with open(filename, 'w') as ofile:
        header = [
            'Variant Accession',
            'Database',
            'Protein Accession',
            'Gene',
            'Location',
            'Amino Acid Change',
            'Amino Acid Change Category',
            'Reference Peptide',
            'Variant Peptide',
            'Strand',
            'Protein Position',
            'Reference Codon',
            'Variant Codon',
            'Variant Nucleotide Position in Codon',
            'Reference Protein Length',
            'Variant Protein Length'
        ]
        ofile.write('\t'.join(header) + '\n')

        t = tqdm(
        variant_peptides,
        desc = 'Writing peptide meta-data file',
        total = len(variant_peptides),
        unit = 'variant peptides'
        )
        for i, vp in enumerate(t):
            if 'reference_residue' in vp.keys() and 'variant_residue' in vp.keys():
                accession = str(vp['variant'])
                protein_accession =vp['reference_transcript'].accession
                var_aa = str(vp['variant_residue'].amino_acid)
                ref_aa = str(vp['reference_residue'].amino_acid)
                amino_acid_change = f'{ref_aa}->{var_aa}'
                amino_acid_category = 'non-synonymous' if ref_aa != var_aa else 'synonymous'
                reference_peptide = get_peptide_string(vp['reference_peptide'], vp['reference_residue'])
                variant_peptide = get_peptide_string(vp['variant_peptide'], vp['variant_residue'])
                reference_codon = get_codon_string(vp['reference_residue'], vp['reference_nucleotide'])
                variant_codon = get_codon_string(vp['variant_residue'], vp['variant_nucleotide'])
                variant_codon_position = get_codon_position(vp['variant_residue'], vp['variant_nucleotide'])
                
                row = [
                    str(vp['variant']),
                    database,
                    vp['reference_transcript'].accession,
                    str(vp['reference_transcript'].gene),
                    'CDS',
                    amino_acid_change,
                    amino_acid_category,
                    reference_peptide,
                    variant_peptide,
                    str(vp['reference_transcript'].strand),
                    str(vp['variant_residue'].position),
                    reference_codon,
                    variant_codon,
                    str(variant_codon_position),
                    str(len(vp['variant_transcript'].protein.sequence)),
                    str(len(vp['reference_transcript'].protein.sequence))
                ]
                ofile.write('\t'.join(row) + '\n')
            elif 'variant_residue' not in vp.keys():
                row = [
                    str(vp['variant']),
                    database,
                    vp['reference_transcript'].accession,
                    str(vp['reference_transcript'].gene),
                    'UTR',
                    'N/A',
                    'N/A',
                    'N/A',
                    'N/A',
                    str(vp['reference_transcript'].strand),
                    'N/A',
                    'N/A',
                    'N/A',
                    'N/A',
                    'N/A',
                    'N/A'
                ]
                ofile.write('\t'.join(row) + '\n')
            else:
                logging.warn(f'variant peptide found without reference peptide\t{i}')



build_peptide_metadata_table(variant_peptides, 'both_isoseq_short', '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/biosurfer/biosurfer/analysis/results/metadata_peptide_both_isoseq_short.tsv')
            