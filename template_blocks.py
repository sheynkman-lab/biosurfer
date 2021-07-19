# bits of code that represent specific functionality from biosurfer


#%%

# extract gencode appris information

import os
from isomodules import isocreate
from isomodules import isofunc

data_dir = '/Users/gloriasheynkman/Documents/research_drive/projects/biosurfer/data'
odir = './results'

if not os.path.exists(odir):
	os.mkdir(odir)

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v38.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v38.pc_transcripts.fa')

# initiate objects
d_gc = isocreate.init_gen_obj_gc(path_gc_gtf)

#TODO - issue is that when reading in objects from gencode, reads in non-protein-coding genes
# write out appris principle transcript names
with open(os.path.join(odir, 'appris_principle_transcript_names.tsv'), 'w') as ofile:
    for symbol, gene in d_gc.items():
        ofile.write('{}\t{}\n'.format(symbol, gene.appris_orf.name))


#%%









#%%