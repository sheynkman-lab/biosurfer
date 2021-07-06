# %%
from IPython.display import display
import os
import pandas as pd
import pickle
from operator import attrgetter
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

genes = (
    'ACSBG2',
    'APOE',
    'ARHGEF1',
    'ARMC6',
    'ATP1A3',
    'CACTIN',
    'CEACAM5',
    'CLASRP',
    'CYP2B6',
    'CYP2S1',
    'DMAC2',
    'ETV2',
    'GCDH',
    'HMG20B',
    'HNRNPUL1',
    'ICAM4',
    'KLF2',
    'KLK3',
    'LSR',
    'MBOAT7',
    'NOSIP',
    'OSCAR',
    'SELENOW',
    'TIMM50',
)
temp_file_path = data_dir + f'gene_dict_{"_".join(sorted(genes))}.p'

try:
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

assert set(genes) <= gd.keys()
aln_grp_dict = dict()

#%%
for gene_name in genes:
    gene = gd[gene_name]
    orfs = sorted(gene, key=lambda orf: (orf is not gene.repr_orf, orf.name))
    aln_grp_dict[gene_name] = isocreatealign.create_and_map_splice_based_align_obj([[orfs[0], orf] for orf in orfs[1:]])

# %%
broken = set()  # FIXME: trying to create Frame object for alignments raises IndexError
force_plotting = False

for gene_name in genes:
    gene = gd[gene_name]
    orfs = sorted(gene, key=lambda orf: (orf is not gene.repr_orf, orf.name))

    fig_path = os.path.join('./data/plots', f'{gene_name}_isoforms.png')
    if force_plotting or not os.path.isfile(fig_path):    
        fig = plt.figure()
        isoplot = isoimage.IsoformPlot(orfs, intron_spacing=30, track_spacing=1.5)
        isoplot.draw()

        try:
            isoplot.draw_frameshifts(aln_grps=aln_grp_dict[gene_name])
        except IndexError:
            broken.add(gene_name)

        fig.set_size_inches(9, 0.5*len(isoplot.orfs))
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300)
        plt.close(fig)
    else:
        aln_grps = isocreatealign.create_and_map_splice_based_align_obj([[orfs[0], orf] for orf in orfs[1:]])

    aln_grp_dict[gene_name] = aln_grps

#%%
# sblocks_dict = {
#     'anchor_orf': [],
#     'other_orf': [],
#     'category': [],
#     'annotation': [],
#     'anchor_first_res': [],
#     'anchor_last_res': [],
#     'other_first_res': [],
#     'other_last_res': []
# }
pblocks_dict = {
    'anchor_orf': [],
    'other_orf': [],
    'category': [],
    'annotation': [],
    'anchor_first_res': [],
    'anchor_last_res': [],
    'other_first_res': [],
    'other_last_res': []
}

for gene_name, aln_grps in aln_grp_dict.items():
    for aln_grp in aln_grps:
        anchor = repr(aln_grp.anchor_orf)
        other = repr(aln_grp.other_orf)
        for i, block in enumerate(aln_grp.alnf.protblocks):
            if block.cat == 'M':
                continue

            annotation = []
            # "first" and "last" refer to exons associated with the pblock itself
            # "prev" and "next" refer to exons associated with the neighboring pblocks, which are type M (match)
            first_alnr, last_alnr = block.first, block.last
            prev_block, prev_exon_anchor, prev_exon_other = None, None, None
            next_block, next_exon_anchor, next_exon_other = None, None, None
            nterm_affected = i == 0
            cterm_affected = i+1 == len(aln_grp.alnf.protblocks)
            if not nterm_affected:
                prev_block = aln_grp.alnf.protblocks[i-1]
                prev_exon_anchor = prev_block.last.res1.exon
                prev_exon_other = prev_block.last.res2.exon
            if not cterm_affected:
                next_block = aln_grp.alnf.protblocks[i+1]
                next_exon_anchor = next_block.first.res1.exon
                next_exon_other = next_block.first.res2.exon

            if block.cat == 'I':
                first_exon = max(first_alnr.res2.exons, key=attrgetter('ord'))
                last_exon = min(last_alnr.res2.exons, key=attrgetter('ord'))

                if not nterm_affected and not cterm_affected:
                    if prev_exon_other is next_exon_other:
                        annotation.append(f'retained intron between exons {prev_exon_anchor.ord} and {next_exon_anchor.ord}')
                    else:
                        number_of_cassette_exons = last_exon.ord - first_exon.ord + 1
                        if prev_exon_other is first_exon:
                            annotation.append(f'exon {prev_exon_anchor.ord} extended by alternative splice donor')
                            number_of_cassette_exons -= 1
                        if next_exon_other is last_exon:
                            annotation.append(f'exon {next_exon_anchor.ord} extended by alternative splice acceptor')
                            number_of_cassette_exons -= 1
                        if number_of_cassette_exons > 0:
                            annotation.append(f'{number_of_cassette_exons} cassette exons inserted between exons {prev_exon_anchor.ord}-{next_exon_anchor.ord}')
                else:
                    if nterm_affected:
                        annotation.append('<N-terminal event>')
                    if cterm_affected:
                        annotation.append('<C-terminal event>')
                    
            elif block.cat == 'D':
                first_exon = max(first_alnr.res1.exons, key=attrgetter('ord'))
                last_exon = min(last_alnr.res1.exons, key=attrgetter('ord'))

                if not nterm_affected and not cterm_affected:
                    if prev_exon_anchor is next_exon_anchor:
                        annotation.append(f'intronized region in exon {prev_exon_anchor.ord}')
                    else:
                        first_cassette_exon_ord = first_exon.ord
                        last_cassette_exon_ord = last_exon.ord
                        if prev_exon_anchor is first_exon:
                            annotation.append(f'exon {prev_exon_anchor.ord} truncated by alternative splice donor')
                            first_cassette_exon_ord += 1
                        if next_exon_anchor is last_exon:
                            annotation.append(f'exon {next_exon_anchor.ord} truncated by alternative splice acceptor')
                            last_cassette_exon_ord -= 1
                        if first_cassette_exon_ord == last_cassette_exon_ord:
                            annotation.append(f'exon {first_cassette_exon_ord} skipped')
                        elif first_cassette_exon_ord < last_cassette_exon_ord:
                            annotation.append(f'exons {first_cassette_exon_ord}-{last_cassette_exon_ord} skipped')
                else:
                    if nterm_affected:
                        annotation.append('<N-terminal event>')
                    if cterm_affected:
                        annotation.append('<C-terminal event>')                
                
            elif block.cat == 'S':
                annotation.append('¯\_(ツ)_/¯')

            pblocks_dict['anchor_orf'].append(anchor)
            pblocks_dict['other_orf'].append(other)
            pblocks_dict['anchor_first_res'].append(first_alnr.res1.idx)
            pblocks_dict['anchor_last_res'].append(last_alnr.res1.idx)
            pblocks_dict['other_first_res'].append(first_alnr.res2.idx)
            pblocks_dict['other_last_res'].append(last_alnr.res2.idx)
            pblocks_dict['category'].append(block.cat)
            pblocks_dict['annotation'].append(', \n'.join(annotation))
        
        # for block in aln_grp.alnf.blocks:
        #     if block.cat == 'M':
        #         continue
        #     anchor_res = (block.first.res1.idx, block.last.res1.idx)
        #     other_res = (block.first.res2.idx, block.last.res2.idx)
        #     sblocks_dict['anchor_orf'].append(anchor)
        #     sblocks_dict['other_orf'].append(other)
        #     sblocks_dict['anchor_first_res'].append(anchor_res[0])
        #     sblocks_dict['anchor_last_res'].append(anchor_res[1])
        #     sblocks_dict['other_first_res'].append(other_res[0])
        #     sblocks_dict['other_last_res'].append(other_res[1])
        #     sblocks_dict['category'].append(block.cat)
        #     sblocks_dict['annotation'].append(None)

pblocks = pd.DataFrame(pblocks_dict)
# sblocks = pd.DataFrame(sblocks_dict)
display(pblocks)
# display(sblocks)
pblocks.to_csv('./data/chr19_annotations.csv', sep='\t', index=False)

# %%
