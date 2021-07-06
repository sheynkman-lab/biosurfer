# %%
from IPython.display import display
import os
import pandas as pd
import pickle
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

#%%
broken = set()  # FIXME: trying to create Frame object for alignments raises IndexError

sblocks_dict = {
    'anchor_orf': [],
    'other_orf': [],
    'category': [],
    'annotation': [],
    'anchor_first_res': [],
    'anchor_last_res': [],
    'other_first_res': [],
    'other_last_res': []
}

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

for gene_name in genes:
    
    gene = gd[gene_name]
    orfs = sorted(gene, key=lambda orf: (orf is not gene.repr_orf, orf.name))

    fig = plt.figure()
    isoplot = isoimage.IsoformPlot(orfs, intron_spacing=30, track_spacing=1.5)
    isoplot.draw()

    try:
        aln_grps = isoplot.draw_frameshifts()
    except IndexError:
        broken.add(gene_name)
        aln_grps = isocreatealign.create_and_map_splice_based_align_obj([[orfs[0], orf] for orf in orfs[1:]])

    fig.set_size_inches(9, 0.5*len(isoplot.orfs))

    for aln_grp in aln_grps:
        anchor = repr(aln_grp.anchor_orf)
        other = repr(aln_grp.other_orf)
        for i, block in enumerate(aln_grp.alnf.protblocks):
            if block.cat == 'M':
                continue

            annotation = []
            first_alnr, last_alnr = block.first, block.last
            prev_block, prev_exon_anchor, prev_exon_other = None, None, None
            next_block, next_exon_anchor, next_exon_other = None, None, None
            
            if block.cat == 'I':
                first_exon = first_alnr.res2.exon
                last_exon = last_alnr.res2.exon
                if i > 0:
                    prev_block = aln_grp.alnf.protblocks[i-1]
                    prev_exon_anchor = prev_block.last.res1.exon
                    prev_exon_other = prev_block.last.res2.exon
                if i+1 < len(aln_grp.alnf.protblocks):
                    next_block = aln_grp.alnf.protblocks[i+1]
                    next_exon_anchor = next_block.first.res1.exon
                    next_exon_other = next_block.first.res2.exon

                prev_exon_is_extended = prev_exon_other is first_exon
                last_exon_is_extended = next_exon_other is last_exon
                if prev_exon_is_extended and last_exon_is_extended:
                    if prev_exon_other is next_exon_other:
                        annotation.append(f'retained intron(s) between exons {prev_exon_anchor.ord} and {next_exon_anchor.ord}')
                    else:
                        pass
                elif prev_exon_is_extended and not last_exon_is_extended:
                    annotation.append(f'exon {prev_exon_anchor.ord} extended by alternate splice donor')
                elif not prev_exon_is_extended and last_exon_is_extended:
                    annotation.append(f'exon {next_exon_anchor.ord} extended by alternate splice acceptor')
                elif prev_block is None:
                    annotation.append('alternative TSS or 5\' UTR')
                elif next_block is None:
                    annotation.append('alternative polyA or 3\' UTR')
                else:
                    annotation.append(f'clean insertion or cassette exons between exons {prev_exon_anchor.ord} and {next_exon_anchor.ord}')
            elif block.cat == 'D':
                first_exon = first_alnr.res1.exon
                last_exon = last_alnr.res1.exon

                first_exon_is_truncated = not first_alnr.res1.is_at_cds_edge
                last_exon_is_truncated = not last_alnr.res1.is_at_cds_edge
                
                if first_exon_is_truncated and last_exon_is_truncated:
                    annotation.append(f'exon {first_exon.ord} truncated')
                    annotation.append(f'exons {first_exon.ord+1}-{last_exon.ord-1} not spliced in')
                    annotation.append(f'exon {last_exon.ord} truncated')
                elif first_exon_is_truncated and not last_exon_is_truncated:
                    annotation.append(f'exon {first_exon.ord} truncated')
                    annotation.append(f'exons {first_exon.ord+1}-{last_exon.ord} not spliced in')
                elif not first_exon_is_truncated and last_exon_is_truncated:
                    annotation.append(f'exons {first_exon.ord}-{last_exon.ord-1} not spliced in')
                    annotation.append(f'exon {last_exon.ord} truncated')
                else:
                    annotation.append(f'exons {first_exon.ord}-{last_exon.ord} not spliced in')
                
            elif block.cat == 'S':
                pass

            pblocks_dict['anchor_orf'].append(anchor)
            pblocks_dict['other_orf'].append(other)
            pblocks_dict['anchor_first_res'].append(first_alnr.res1.idx)
            pblocks_dict['anchor_last_res'].append(last_alnr.res1.idx)
            pblocks_dict['other_first_res'].append(first_alnr.res2.idx)
            pblocks_dict['other_last_res'].append(last_alnr.res2.idx)
            pblocks_dict['category'].append(block.cat)
            pblocks_dict['annotation'].append(', '.join(annotation))
        
        for block in aln_grp.alnf.blocks:
            if block.cat == 'M':
                continue
            anchor_res = (block.first.res1.idx, block.last.res1.idx)
            other_res = (block.first.res2.idx, block.last.res2.idx)
            sblocks_dict['anchor_orf'].append(anchor)
            sblocks_dict['other_orf'].append(other)
            sblocks_dict['anchor_first_res'].append(anchor_res[0])
            sblocks_dict['anchor_last_res'].append(anchor_res[1])
            sblocks_dict['other_first_res'].append(other_res[0])
            sblocks_dict['other_last_res'].append(other_res[1])
            sblocks_dict['category'].append(block.cat)
            sblocks_dict['annotation'].append(None)

    plt.show()

#%%

pblocks = pd.DataFrame(pblocks_dict)
sblocks = pd.DataFrame(sblocks_dict)
display(pblocks)
# display(sblocks)

# # %%
# for gene_name in broken:
#     gene = gd[gene_name]
#     aln_grps = isocreatealign.create_and_map_splice_based_align_obj([[gene.repr_orf, orf] for orf in gene.other_orfs])
#     for aln_grp in aln_grps:
#         isocreatefeat.create_and_map_frame_objects(aln_grp)
# %%
