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

data_dir = "/home/redox/sheynkman-lab/gencode"

# filepaths
path_chr_gtf = os.path.join(data_dir, "chr19.gtf")
path_chr_fa = os.path.join(data_dir, "chr19.fa")
path_transcripts_fa = os.path.join(data_dir, "gencode.v38.pc_transcripts.fa")
path_hg38_fa = os.path.join(data_dir, 'GRCh38.primary_assembly.genome.fa')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_transcripts_fa)
chr_dict = isofunc.load_hg38(path_chr_fa)
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
genes = ['HMG20B']
gd = isocreate.init_gen_obj_gc(path_chr_gtf, gene_list=genes)
gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs)
gd = isocreate.create_and_link_junct_and_ss_objs(gd, chr_dict)
# gd = isocreate.create_and_map_domains(gd, domains)

# %%
reload(isoimage)
isoimage.render_iso_image(list(gd[genes[0]].orfs))
