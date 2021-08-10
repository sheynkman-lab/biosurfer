#%%
import traceback
from operator import attrgetter

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
        aln_dict[gene] = [TranscriptBasedAlignment(isoform_list[0], other) for other in isoform_list[1:]]
    except Exception as e:
        print(f'----------------\n{gene}')
        traceback.print_exc()

#%%
# example of frameshift on plus strand (f-category) and complex split codon alignment (x-category)
aln = TranscriptBasedAlignment(proteins['TANGO2-201'], proteins['TANGO2-207'])
print(repr(aln))
print(aln.full)

#%%
# example of frameshift on minus strand (f-category)
aln = TranscriptBasedAlignment(proteins['RBFOX2-201'], proteins['RBFOX2-202'])
print(repr(aln))
print(aln.full)

#%%
# example of "edge mismatch" (e-category)
aln = TranscriptBasedAlignment(proteins['MAPK12-201'], proteins['MAPK12-202'])
print(repr(aln))
print(aln.full)

#%%
# example of single-nt overlap between two codons (treated as separate alignment pairs)
aln = TranscriptBasedAlignment(proteins['BID-201'], proteins['BID-202'])
print(repr(aln))
print(aln.full)

# %%
