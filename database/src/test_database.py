#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from alignments import TranscriptBasedAlignment
from sqlalchemy import select
# from itertools import combinations, islice

def get_gene_protein_isoforms(gene_name):
    gene = db_session.execute(select(Gene).filter(Gene.name == gene_name)).one()[Gene]
    return {transcript.name: transcript.orfs[0].protein for transcript in gene.transcripts}

#%%
proteins = get_gene_protein_isoforms('TANGO2')
proteins.update(get_gene_protein_isoforms('RBFOX2'))
proteins.update(get_gene_protein_isoforms('MAPK12'))
proteins.update(get_gene_protein_isoforms('BID'))

#%%
# example of frameshift on plus strand (f-category)
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
