#%%
import os
import traceback
from operator import attrgetter
from warnings import filterwarnings

import matplotlib.pyplot as plt
import pandas as pd
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning
from sqlalchemy import select

from alignments import TranscriptBasedAlignment
from database import db_session
from models import Transcript, Exon, Gene, Protein, Transcript
from plotting import IsoformPlot

filterwarnings("ignore", category=MatplotlibDeprecationWarning)


def get_transcripts_sorted_by_appris(gene):
    return sorted(gene.transcripts, key=attrgetter('appris'))


def get_protein_isoforms(gene):
    return {transcript.name: transcript.orfs[0].protein for transcript in get_transcripts_sorted_by_appris(gene) if transcript.orfs}


#%%
gene_list = (
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
        transcript_list = get_transcripts_sorted_by_appris(gene)
        transcripts[name] = transcript_list
        isoforms = {transcript.name: transcript.orfs[0].protein for transcript in transcript_list if transcript.orfs}
        proteins.update(isoforms)
        isoform_list = list(isoforms.values())
        anchor = isoform_list[0]
        for other in isoform_list[1:]:
            aln_dict[other.orf.transcript.name] = TranscriptBasedAlignment(anchor, other)

# %%
noncoding_transcripts = {transcript for gene in genes.values() for transcript in gene.transcripts if not transcript.orfs}

# %%
broken = set()
force_plotting = True

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
            print(e)
        else:
            fig.set_size_inches(9, 0.5*len(isoplot.transcripts))
            plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
            print('saved '+fig_path)
        finally:
            plt.close(fig)

# %%
for aln in aln_dict.values():
    aln.annotate()

all_annotations = pd.DataFrame.from_records(
    ((str(aln.anchor), str(aln.other), str(pblock.category), pblock.annotation) for aln in aln_dict.values() for pblock in aln.protein_blocks if pblock.annotation),
    columns=('anchor', 'other', 'category', 'annotation')
)
display(all_annotations)
all_annotations.to_csv('chr22_annotations.tsv', sep='\t')
