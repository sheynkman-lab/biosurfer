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

from biosurfer.core.alignments import TranscriptBasedAlignment
from biosurfer.core.models import ORF, Gene, Protein, db_session
from biosurfer.core.models import GencodeTranscript as Transcript
from biosurfer.plots.plotting import IsoformPlot

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

#%%
# proteins = db_session.query(Protein).join(ORF, Transcript).where(Transcript.basic).order_by(func.random()).limit(20).all()
# genes = {protein.gene for protein in proteins}

gene_names = (
    'ANKRD6',
    'APOBEC3B',
    'BID',
    'CAPN3',
    'CCDC178',
    'CSTF3',
    'DERL3',
    'ETV2',
    'HNF4A',
    'HOMEZ',
    'MECP2',
    'MEGF10',
    'NAXD',
    'POLR2F',
    'RAD51AP1',
    'SLC2A11',
    'TANGO2',
    'TMPO',
    'TTC23L',
    'WIPF1',
)
genes = [Gene.from_name(name) for name in gene_names]
display(genes)

#%%
force_plotting = False
aln_dict = dict()
for gene in genes:
    transcript_list = sorted((tx for tx in gene.transcripts if tx.basic), key=attrgetter('appris'))
    isoforms = [transcript.orfs[0].protein for transcript in transcript_list]
    anchor = isoforms[0]
    for other in isoforms[1:]:
        aln_dict[other.orf.transcript.name] = TranscriptBasedAlignment(anchor, other)

    fig_path = f'../data/plots/{gene.name}_isoforms.png'
    if force_plotting or not os.path.isfile(fig_path):
        fig = plt.figure()
        isoplot = IsoformPlot(transcript_list)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        k = 1
        while 2*k*max(orf.stop - orf.start for tx in transcript_list for orf in tx.orfs) < max(tx.stop - tx.start for tx in transcript_list):
            k *= 2
        fig.set_size_inches(9*k, 0.5*len(isoplot.transcripts))
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
        print('saved '+fig_path)
        plt.close(fig)

# %%
annotations_output = 'sample_annotations.tsv'
with open(annotations_output, 'w') as f:
    writer = csv.writer(f, delimiter='\t', quotechar='"')
    writer.writerow(['anchor', 'other', 'category', 'region', 'event', 'annotation'])
    for aln in aln_dict.values():
        try:
            aln.annotate()
        except Exception as e:
            print(aln)
            traceback.print_exc()
        for pblock in aln.protein_blocks:
            if pblock.annotation:
                writer.writerow([str(aln.anchor), str(aln.other), str(pblock.category), str(pblock.region), pblock.event, pblock.annotation])

# %%
alignments_output = 'sample_alignments.txt'
all_full = '\n\n\n'.join(str(aln)+'\n'+aln.full for aln in aln_dict.values())
with open(alignments_output, 'w') as f:
    f.write(all_full)

# %%
