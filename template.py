# starter template to create iso-objects and "surf" analysis

#%%

# during development, need to reload modules
from importlib import reload

import os
from isomodules import isocreate
from isomodules import isocreatealign
from isomodules import isocreatefeat
from isomodules import isoclass
from isomodules import isofunc
from isomodules import isofeature
from isomodules import isomap
from isomodules import isoimage
from isomodules import isogroup
from isomodules import isoalign
from isomodules import isowrangle
from isomodules import isocompare
from bx.intervals.intersection import Interval, IntervalTree

data_dir = './data'
odir = 'results'

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v38.annotation.gtf.toy')
path_gc_fa = os.path.join(data_dir, 'gencode.v38.pc_transcripts.fa.toy')
path_gc_domain_toy = os.path.join(data_dir, 'gencode.v30.domains.toy')
path_pfam_gc = os.path.join(data_dir, 'pfam/2019-07-04_HMMER_domain_mappings_to_GS_fasta_file.txt')
path_pfam_names = os.path.join(data_dir, 'pfam/pfam_a_names.tsv')
path_hg38_fa = os.path.join(data_dir, 'GRCh38_canon.fa')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
hg38_dict = isofunc.load_hg38(path_hg38_fa)
domain_dict = isofunc.load_domain_mappings(path_gc_domain_toy, path_pfam_names)

def make_directory(dname):
    if not os.path.exists(dname):
        os.mkdir(dname)

make_directory('output_isoimages')

genes = ['NFYA', 'PAX5']

d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs)
d_gc_w_seqs_w_juncs = isocreate.create_and_link_junct_and_ss_objs(d_gc_w_seqs, hg38_dict)

# map domains and write out text images 
all_domains = []
oname = os.path.join(odir, 'domain_mappings.txt')
with open(oname, 'w') as ofile:
    for gene_name, gene in d_gc_w_seqs_w_juncs.items():
        for orf in gene:
            domain_infos = domain_dict[orf.name]
            for domain_info in domain_infos: 
                domain = isocreatefeat.create_and_map_domain(orf, domain_info)
                all_domains.append(domain)
                ofile.write('\n')
                ofile.write(domain.full)

print(all_domains)
gd = d_gc_w_seqs_w_juncs


isoimage.render_iso_image(gd['NFYA'].orfs)
