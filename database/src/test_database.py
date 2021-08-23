#%%
import csv
import os
import traceback
from operator import attrgetter
from warnings import filterwarnings

import matplotlib.pyplot as plt
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

from alignments import TranscriptBasedAlignment
from models import Exon, Gene, Protein, Transcript
from plotting import IsoformPlot

filterwarnings("ignore", category=MatplotlibDeprecationWarning)


def get_transcripts_sorted_by_appris(gene):
    return sorted(gene.transcripts, key=attrgetter('appris'))


def get_protein_isoforms(gene):
    return {transcript.name: transcript.orfs[0].protein for transcript in get_transcripts_sorted_by_appris(gene) if transcript.orfs}


#%%
chr22_genes = (
    'APOBEC3B',
    'BID',
    'CABIN1',
    'CHEK2',
    'DERL3',
    'EWSR1',
    'GGT1',
    'GUCD1',
    'INPP5J',
    'LARGE1',
    'MAPK12',
    'MICAL3',
    'NF2',
    'PISD',
    'POLR2F',
    'RAC2',
    'RBFOX2',
    'SEPTIN5',
    'SERHL2',
    'SHANK3',
    'SLC2A11',
    'SMTN',
    'SPECC1L',
    'SYNGR1',
    'TANGO2',
)
chr19_genes = (
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
    # 'SELENOW',
    'TIMM50',
)
gene_list = chr19_genes
annotations_output = 'chr19_annotations.tsv'
alignments_output = 'chr19_alignments.txt'
genes = dict()
transcripts = dict()
proteins = dict()
aln_dict = dict()

for name in gene_list:
    try:
        gene = Gene.from_name(name)
    except Exception as e:
        print(f'----------------\ncould not get {name}')
        traceback.print_exc()
    else:
        genes[name] = gene
        transcript_list = [tx for tx in get_transcripts_sorted_by_appris(gene) if tx.basic and tx.orfs]
        transcripts[name] = transcript_list
        try:
            isoforms = {transcript.name: transcript.orfs[0].protein for transcript in transcript_list}
            proteins.update(isoforms)
            isoform_list = list(isoforms.values())
            anchor = isoform_list[0]
            for other in isoform_list[1:]:
                aln_dict[other.orf.transcript.name] = TranscriptBasedAlignment(anchor, other)
        except Exception as e:
            print(f'----------------\ncould not get proteins for {name}: {e}')

# %%
noncoding_transcripts = {transcript for gene in genes.values() for transcript in gene.transcripts if not transcript.orfs}

# %%
broken = set()
force_plotting = False

for name, tx_list in transcripts.items():
    fig_path = f'../../data/plots/{name}_isoforms.png'
    if force_plotting or not os.path.isfile(fig_path):
        try:
            fig = plt.figure()
            isoplot = IsoformPlot(tx_list)
            isoplot.draw_all_isoforms()
            isoplot.draw_frameshifts()
        except Exception as e:
            broken.add(name)
            print(f'----------------\ncould not plot {name}')
            traceback.print_exc()
        else:
            fig.set_size_inches(9, 0.5*len(isoplot.transcripts))
            plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
            print('saved '+fig_path)
        finally:
            plt.close(fig)

# %%
with open(annotations_output, 'w') as f:
    writer = csv.writer(f, delimiter='\t', quotechar='"')
    writer.writerow(['anchor', 'other', 'category', 'region', 'annotation'])
    for aln in aln_dict.values():
        aln.annotate()
        for pblock in aln.protein_blocks:
            if pblock.annotation:
                writer.writerow([str(aln.anchor), str(aln.other), str(pblock.category), str(pblock.region), pblock.annotation])

# %%
all_full = '\n\n\n'.join(str(aln)+'\n'+aln.full for aln in aln_dict.values())
with open(alignments_output, 'w') as f:
    f.write(all_full)

# %%
