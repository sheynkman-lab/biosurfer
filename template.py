# title           :driver.py
# description     :Starting template to iso-objs and do analysis.
# author          :Gloria Sheynkman
# ==============================================================================
# test

# %%
import os
from imp import reload
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

data_dir = "/home/redox/sheynkman-lab/gencode"

# filepaths
path_chr22_gtf = os.path.join(data_dir, "chr22.gtf")
path_transcripts_fa = os.path.join(data_dir, "gencode.v38.pc_transcripts.fa")
path_hg38_fa = os.path.join(data_dir, 'GRCh38.primary_assembly.genome.fa')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_transcripts_fa)
# hg38_dict = isofunc.load_hg38(path_hg38_fa)
# domains = isofunc.load_domain_mappings(path_pfam_gc, path_pfam_names)
# aln_data = isofunc.load_clustal_alignments(path_gc_aln)
# aln_blocks = isofunc.load_alignment_block_calls(aln_data, path_gc_aln_blocks)
# tf_sachi, tf_gloria, tf_lambert = isofunc.load_tf_list(path_tf_list_sachi,
#                                      path_tf_list_gloria, path_tf_list_lambert)
# tfs = isofunc.load_man_correct_gc30_tf_genenames(path_gc_corr_tf_list)
# appris_orfs = isofunc.load_appris_principle_isonames(path_appris)
# isoacc_map = isofunc.load_6k_isoacc_map(path_isoacc)

# %%

def make_directory(dname):
    if not os.path.exists(dname):
        os.mkdir(dname)

make_directory('output_isoimages')


reload(isocreate)
reload(isoclass)
reload(isofunc)
reload(isofeature)
reload(isomap)
reload(isoimage)
reload(isogroup)
reload(isoalign)
reload(isowrangle)

genes = ['NFYA', 'PAX5']

d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs)
d_gc_w_seqs_w_juncs = isocreate.create_and_link_junction_related_obj(d_gc_w_seqs, hg38_dict)
d_gc_w_seqs_w_juncs_w_doms = isocreate.create_and_map_domains(d_gc_w_seqs_w_juncs, domains)
gd = d_gc_w_seqs_w_juncs_w_doms  # gen_obj dict

# %%
genes = ['MAPK1', 'MAPK12']
d_gc = isocreate.init_gen_obj_gc(path_chr22_gtf, gene_list=genes)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs)

# %%

reload(isoimage)
isoimage.render_iso_image(gd['NFYA'].orfs)


genes = ['NFYA', 'PAX5']
gd = isocreate.init_gen_obj(path_6k_gtf, genes)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs_6k)





# %%

# clustal alignments into OOP objects
grps = isocreate.create_and_map_alignments(gd, aln_data, aln_blocks)
for grp in grps:
    print(grp.alnf)
    grp.enter()
    print(grp.anchor_orf.current_grp)
    alnf = grp.anchor_orf.aln
    print(''.join([alnr.alnb.cat[1] if alnr not in [alnr.alnb.first, alnr.alnb.last] else '|' for alnr in alnf.chain]))
    print(''.join([alnr.alnsb.cat[1] if alnr not in [alnr.alnsb.first, alnr.alnsb.last] else '|' for alnr in alnf.chain]))
    print(''.join([str(alnr.alnsb.cds1.ord)[0] for alnr in alnf.chain]))
    print(''.join([alnr.res1.aa for alnr in alnf.chain]))
    print(''.join([alnr.cat for alnr in alnf.chain]))
    print(''.join([alnr.res2.aa for alnr in alnf.chain]))
    print(''.join([str(alnr.alnsb.cds2.ord)[0] for alnr in alnf.chain]))
    print('\n')
    grp.exit()
