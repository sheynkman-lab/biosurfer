from time import time

from biosurfer.core.populate_database import (FastaHeaderFields,
                                              load_gencode_gtf,
                                              load_pacbio_gtf,
                                              load_sqanti_classifications,
                                              load_transcript_fasta,
                                              load_translation_fasta)

path = '../data/'
gencode_gtf = 'gencode/gencode.v38.annotation.gtf'
gencode_tx = 'gencode/gencode.v38.transcripts.fa'
gencode_tl = 'gencode/gencode.v38.pc_translations.fa'
pacbio_gtf = 'bone/bone_cds_high_confidence.gtf'
pacbio_tx = 'bone/filtered_bone_corrected.fasta'
pacbio_tl = 'bone/bone_hybrid.fasta'
sqanti = 'bone/SQANTI3_results_attempt2_classification.tsv'

def get_ids_from_gencode_fasta(header: str):
    fields = [field for field in header.split('|') if field and not field.startswith(('UTR', 'CDS'))]
    transcript_id = next((field for field in fields if field.startswith('ENST')))
    protein_id = next((field for field in fields if field.startswith('ENSP')), None)
    return FastaHeaderFields(transcript_id, protein_id)

def get_ids_from_pacbio_fasta(header: str):
    return FastaHeaderFields(header, None)

def get_ids_from_hybrid_fasta(header: str):
    fields = header.split('|')
    return FastaHeaderFields(fields[1], fields[1] + ':PROT1')

def skip_par_y(header: str):  # these have duplicate ENSEMBL accessions and that makes SQLAlchemy very sad
    return 'PAR_Y' in header

def skip_gencode(header: str):
    return header.startswith('gc')

#%%
start = time()
load_gencode_gtf(path + gencode_gtf, overwrite=False)
load_transcript_fasta(path + gencode_tx, get_ids_from_gencode_fasta, skip_par_y)
load_translation_fasta(path + gencode_tl, get_ids_from_gencode_fasta, skip_par_y)

#%%
load_pacbio_gtf(path + pacbio_gtf, overwrite=False)
load_transcript_fasta(path + pacbio_tx, get_ids_from_pacbio_fasta)
load_translation_fasta(path + pacbio_tl, get_ids_from_hybrid_fasta, skip_gencode)

# %%
load_sqanti_classifications(path + sqanti)

# %%
end = time()
print(f'Total time: {end - start:.3g} s')
