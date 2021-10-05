from time import time

from biosurfer.core.database import Database, DB_GENCODE
from biosurfer.core.helpers import get_ids_from_gencode_fasta, skip_par_y

path = '../data/'
gencode_gtf = 'gencode/gencode.v38.annotation.gtf'
gencode_tx = 'gencode/gencode.v38.transcripts.fa'
gencode_tl = 'gencode/gencode.v38.pc_translations.fa'
gencode_doms = 'gencode/grch38-pfam-mappings.tsv'
pfam_dom_names = 'gencode/pfam_a_names.tsv'

#%%
db = Database(DB_GENCODE)

start = time()
db.load_gencode_gtf(path + gencode_gtf, overwrite=True)
db.load_transcript_fasta(path + gencode_tx, get_ids_from_gencode_fasta, skip_par_y)
db.load_translation_fasta(path + gencode_tl, get_ids_from_gencode_fasta, skip_par_y)
db.load_domain_mappings(path + gencode_doms, path + pfam_dom_names)
end = time()
print(f'Total time: {end - start:.3g} s')
