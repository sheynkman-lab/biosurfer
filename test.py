# %%
from IPython.display import display
import os
import pandas as pd
from importlib import reload
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

import matplotlib.pyplot as plt
from matplotlib._api.deprecation import MatplotlibDeprecationWarning
from warnings import filterwarnings
filterwarnings("ignore", category=MatplotlibDeprecationWarning)

import pickle

# data_dir = '/home/redox/sheynkman-lab/gencode'
data_dir = './data/biosurfer_demo_data/'

# filepaths
# jared's fpaths
path_chr_gtf = os.path.join(data_dir, 'chr19.gtf')
path_chr_fa = os.path.join(data_dir, 'chr19.fa')
path_transcripts_fa = os.path.join(data_dir, 'gencode.v38.pc_transcripts.fa')

# gloria's fpaths
# path_chr_gtf = os.path.join(data_dir, 'chr19.gtf')
# path_chr_fa = os.path.join(data_dir, 'chr19.fa')
# path_transcripts_fa = os.path.join(data_dir, 'gencode.v38.pc_transcripts.fa')
# path_hg38_fa = os.path.join(data_dir, 'GRCh38.primary_assembly.genome.fa')

genes = ('HMG20B', 'DMAC2', 'APOE', 'KLF2', 'TIMM50', 'HRC', 'GCDH')
temp_file_path = f'data/gene_dict_{"_".join(sorted(genes))}.p'

try:
    # raise IOError()
    with open(temp_file_path, 'rb') as f:
        gd = pickle.load(f)
        print("loaded dict from pickle file")
except IOError:
    print("loading dict from gtf and fasta files")
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

    gd = isocreate.init_gen_obj_gc(path_chr_gtf, gene_list=genes)
    gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs)
    gd = isocreate.create_and_link_junct_and_ss_objs(gd, chr_dict)
    # gd = isocreate.create_and_map_domains(gd, domains)

    with open(temp_file_path, 'wb') as f:
        pickle.dump(gd, f)

#%%

# for junc in goi.repr_orf.juncs:
#     print(junc.up_exon, junc, junc.dn_exon)
#     if junc.len > 1000:
#         print(junc.up_exon.first.res)

#%%
alt_regs_dict = {
    'anchor': [],
    'other': [],
    'anchor_first': [],
    'anchor_last': [],
    'other_first': [],
    'other_last': [],
    'category': []
}

broken = ('GCDH', 'TIMM50')  # FIXME: trying to create Frame object for alignments raises IndexError
for gene_name in sorted(genes):
    if gene_name in broken:
        continue
    gene = gd[gene_name]
    orfs = sorted(gene, key=lambda orf: (orf is not gene.repr_orf, orf.name))

    fig = plt.figure()
    isoplot = isoimage.IsoformPlot(orfs, intron_spacing=30, track_spacing=1.5)
    isoplot.draw()
    aln_grps = isoplot.draw_frameshifts()

    # print('anchor      anchor res     other       other res      category ')
    for aln_grp in aln_grps:
        anchor = repr(aln_grp.anchor_orf)
        other = repr(aln_grp.other_orf)
        for block in aln_grp.alnf.blocks:
            if block.cat == 'M':
                continue
            anchor_res = (block.first.res1.idx, block.last.res1.idx)
            other_res = (block.first.res2.idx, block.last.res2.idx)
            # print(f'{anchor:12}{str(anchor_res):15}{other:12}{str(other_res):15}{block.cat}')
            alt_regs_dict['anchor'].append(anchor)
            alt_regs_dict['other'].append(other)
            alt_regs_dict['anchor_first'].append(anchor_res[0])
            alt_regs_dict['anchor_last'].append(anchor_res[1])
            alt_regs_dict['other_first'].append(other_res[0])
            alt_regs_dict['other_last'].append(other_res[1])
            alt_regs_dict['category'].append(block.cat)

    plt.show()

alt_regs = pd.DataFrame(alt_regs_dict)
display(alt_regs)



# %%
