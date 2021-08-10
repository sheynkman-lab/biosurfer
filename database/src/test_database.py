#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from alignments import TranscriptBasedAlignment
from sqlalchemy import select
# from itertools import combinations

def get_gene_protein_isoforms(gene_name):
    gene = db_session.execute(select(Gene).filter(Gene.name == gene_name)).one()[Gene]
    return {transcript.name: transcript.orfs[0].protein for transcript in gene.transcripts if transcript.orfs}

#%%
genes = (
    'APOBEC3B',
    'BID',
    'CHEK2',
    'EWSR1',
    'LARGE1',
    'MAPK12',
    'MICAL3',
    'MRTFA',
    'PISD',
    'RAC2',
    'RBFOX2',
    'SEPTIN5',
    'SHANK3',
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
        print(f'{gene}: {repr(e)}')

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
