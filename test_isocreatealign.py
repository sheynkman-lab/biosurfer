# %%
import os
import pickle
from warnings import filterwarnings

from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

from isomodules import (helpers, isoalign, isoclass, isocreate, isocreatealign,
                        isocreatefeat, isofeature, isogroup, isoimage)

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

def get_gene_dictionary(gene_list):
    temp_file_path = data_dir + f'gene_dict_{"_".join(sorted(gene_list))}.p'
    try:
        with open(temp_file_path, 'rb') as f:
            gd = pickle.load(f)
            print("loaded dict from pickle file")
    except IOError:
        print("loading dict from gtf and fasta files")
        # load data
        orf_seqs = helpers.gc_fasta_to_orf_seq_dict(path_transcripts_fa)
        chr_dict = helpers.load_hg38(path_chr_fa)
        # hg38_dict = helpers.load_hg38(path_hg38_fa)
        # domains = helpers.load_domain_mappings(path_pfam_gc, path_pfam_names)
        # aln_data = helpers.load_clustal_alignments(path_gc_aln)
        # aln_blocks = helpers.load_alignment_block_calls(aln_data, path_gc_aln_blocks)
        # tf_sachi, tf_gloria, tf_lambert = helpers.load_tf_list(path_tf_list_sachi,
        #                                      path_tf_list_gloria, path_tf_list_lambert)
        # tfs = helpers.load_man_correct_gc30_tf_genenames(path_gc_corr_tf_list)
        # appris_orfs = helpers.load_appris_principle_isonames(path_appris)
        # isoacc_map = helpers.load_6k_isoacc_map(path_isoacc)

        gd = isocreate.init_gen_obj_gc(path_chr_gtf, gene_list=gene_list)
        gd = isocreate.create_and_link_seq_related_obj(gd, orf_seqs)
        gd = isocreate.create_and_link_junct_and_ss_objs(gd, chr_dict)
        # gd = isocreate.create_and_map_domains(gd, domains)

        with open(temp_file_path, 'wb') as f:
            pickle.dump(gd, f)

    assert set(gene_list) <= gd.keys()
    return gd

gd = get_gene_dictionary(genes)
aln_grp_dict = dict()

# %%
def is_complete(orf):
    start_nf = hasattr(orf, 'cds_start_nf') and orf.cds_start_nf
    end_nf = hasattr(orf, 'cds_end_nf') and orf.cds_end_nf
    return not start_nf and not end_nf

print('filtering out incomplete isoforms')
complete_orfs = {gene_name: sorted(filter(is_complete, gene), key=lambda orf: (orf is not gene.repr_orf, orf.name)) for gene_name, gene in sorted(gd.items())}
print('creating alignment objects')
for gene_name, orfs in complete_orfs.items():
    print(f'\t{gene_name}')
    # aln_grp_dict[gene_name] = isocreatealign.create_and_map_splice_based_align_obj([[orfs[0], orf] for orf in orfs[1:]])
    aln_grp_dict[gene_name] = [isocreatealign.get_splice_aware_isoform_alignment(orfs[0], other_orf) for other_orf in orfs[1:]]

# %%
with open('data/sample_alignment_repr_2.txt', 'w') as f:
    f.write('\n'.join(aln_grp.alnf.full for aln_grps in aln_grp_dict.values() for aln_grp in aln_grps))
