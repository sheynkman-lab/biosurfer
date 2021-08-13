#%%
import traceback
from operator import attrgetter

from IPython.display import display
from sqlalchemy import select

from alignments import TranscriptBasedAlignment
from database import db_session
from models import ORF, Exon, Gene, Protein, Transcript

def get_gene_protein_isoforms(gene_name):
    gene = db_session.execute(select(Gene).filter(Gene.name == gene_name)).one()[Gene]
    return {transcript.name: transcript.orfs[0].protein for transcript in sorted(gene.transcripts, key=attrgetter('appris')) if transcript.orfs}

#%%
genes = (
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
    'MRTFA',
    'NF2',
    'PISD',
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
proteins = dict()
aln_dict = dict()
for gene in genes:
    try:
        isoforms = get_gene_protein_isoforms(gene)
        proteins.update(isoforms)
        isoform_list = list(isoforms.values())
        anchor = isoform_list[0]
        for other in isoform_list[1:]:
            aln_dict[(anchor.orf.transcript.name, other.orf.transcript.name)] = TranscriptBasedAlignment(anchor, other)
    except Exception as e:
        print(f'----------------\n{gene}')
        traceback.print_exc()

# %%
alns = (
    ('TANGO2-201', 'TANGO2-207'),
    ('MAPK12-201', 'MAPK12-202'),
    ('BID-208', 'BID-202'),
)

for anchor, other in alns:
    aln = aln_dict[anchor, other]
    aln.annotate()
    print(f'{aln}\n{aln.full}')
    display(aln.transcript_blocks)
    display(aln.protein_blocks)
    for pblock in aln.protein_blocks:
        print(pblock.annotation)

# %%
