#%%
import csv
import os
import traceback
from operator import attrgetter
from warnings import filterwarnings

import matplotlib.pyplot as plt
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning
from sqlalchemy.sql.expression import and_, func

from alignments import TranscriptBasedAlignment
from models import ORF, Gene, Protein, db_session
from models import GencodeTranscript as Transcript
from plotting import IsoformPlot

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

#%%
# proteins = db_session.query(Protein).join(ORF, Transcript).where(Transcript.basic).order_by(func.random()).limit(20).all()
# genes = {protein.gene for protein in proteins}

gene_names = (
    'ANKRD6',
    'CAPN3',
    'CCDC178',
    'CSTF3',
    'H1-0',
    'HNF4A',
    'HOMEZ',
    'LZIC',
    'MECP2',
    'MEGF10',
    'NAXD',
    'RAD51AP1',
    'RALGPS2',
    'S100A13',
    'TMPO',
    'TSPAN12',
    'TTC23L',
    'VPS53',
    'WIPF1',
    'ZC3H14'
)
genes = [Gene.from_name(name) for name in gene_names]
display(genes)

#%%
annotations_output = 'sample_annotations.tsv'
alignments_output = 'sample_alignments.txt'
aln_dict = dict()
force_plotting = False

for gene in genes:
    transcript_list = sorted((tx for tx in gene.transcripts if tx.basic), key=attrgetter('appris'))
    isoforms = [transcript.orfs[0].protein for transcript in transcript_list]
    anchor = isoforms[0]
    for other in isoforms[1:]:
        aln_dict[other.orf.transcript.name] = TranscriptBasedAlignment(anchor, other)

    fig_path = f'../../data/plots/{gene.name}_isoforms.png'
    if force_plotting or not os.path.isfile(fig_path):
        fig = plt.figure()
        isoplot = IsoformPlot(transcript_list)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        fig.set_size_inches(9, 0.5*len(isoplot.transcripts))
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
        print('saved '+fig_path)
        plt.close(fig)

# %%
with open(annotations_output, 'w') as f:
    writer = csv.writer(f, delimiter='\t', quotechar='"')
    writer.writerow(['anchor', 'other', 'category', 'region', 'event', 'annotation'])
    for aln in aln_dict.values():
        aln.annotate()
        for pblock in aln.protein_blocks:
            if pblock.annotation:
                writer.writerow([str(aln.anchor), str(aln.other), str(pblock.category), str(pblock.region), pblock.event, pblock.annotation])

# %%
all_full = '\n\n\n'.join(str(aln)+'\n'+aln.full for aln in aln_dict.values())
with open(alignments_output, 'w') as f:
    f.write(all_full)

# %%
