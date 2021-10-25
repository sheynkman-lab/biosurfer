from time import time

from biosurfer.core.database import Database
from biosurfer.core.helpers import get_ids_from_gencode_fasta, skip_par_y, FastaHeaderFields

path = '../data/'
gencode_gtf = 'gencode/gencode.v38.basic.annotation.gtf'
gencode_tx = 'gencode/gencode.v38.pc_transcripts.fa'
gencode_tl = 'gencode/gencode.v38.pc_translations.fa'
gencode_doms = 'gencode/grch38-pfam-mappings.tsv'
pfam_dom_info = 'gencode/pfamA.tsv'
pacbio_gtf = 'bone/bone_with_cds.gtf'
pacbio_tx = 'bone/filtered_bone_corrected.fasta'
pacbio_tl = 'bone/bone_orf_refined.fasta'
sqanti = 'bone/SQANTI3_results_attempt2_classification.tsv'

def get_ids_from_pacbio_fasta(header: str):
    return FastaHeaderFields(header, None)

def get_ids_from_hybrid_fasta(header: str):
    fields = header.split('|')
    return FastaHeaderFields(fields[1], fields[1] + ':PROT1')

def skip_gencode(header: str):
    return header.startswith('gc')

#%%
db = Database('bone')
start = time()
db.load_gencode_gtf(path + gencode_gtf, overwrite=False)
db.load_transcript_fasta(path + gencode_tx, get_ids_from_gencode_fasta, skip_par_y)
db.load_translation_fasta(path + gencode_tl, get_ids_from_gencode_fasta, skip_par_y)
db.load_domains(path + pfam_dom_info, overwrite=True)
db.load_domain_mappings(path + gencode_doms, overwrite=True)
db.load_pacbio_gtf(path + pacbio_gtf)
db.load_transcript_fasta(path + pacbio_tx, get_ids_from_pacbio_fasta)
db.load_translation_fasta(path + pacbio_tl, get_ids_from_hybrid_fasta, skip_gencode)
db.load_sqanti_classifications(path + sqanti)
end = time()
print(f'Total time: {end - start:.3g} s')

