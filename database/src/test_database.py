#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from alignments import ProteinAlignment
from alignments import TranscriptLevelEvent as TLE
from sqlalchemy import select
from itertools import combinations, islice

def get_gene_protein_isoforms(gene_name):
    gene = db_session.execute(select(Gene).filter(Gene.name == gene_name)).one()[Gene]
    return {transcript.name: transcript.orfs[0].protein for transcript in gene.transcripts}

#%%
proteins = get_gene_protein_isoforms('TANGO2')


aln = ProteinAlignment(proteins['TANGO2-201'], proteins['TANGO2-202'])
print(repr(aln))
print(aln.full)

#%%
proteins = get_gene_protein_isoforms('RBFOX2')

aln = ProteinAlignment(proteins['RBFOX2-201'], proteins['RBFOX2-202'])
print(repr(aln))
print(aln.full)

#%%
proteins = get_gene_protein_isoforms('MAPK12')

aln = ProteinAlignment(proteins['MAPK12-201'], proteins['MAPK12-202'])
print(repr(aln))
print(aln.full)

#%%
proteins = get_gene_protein_isoforms('BID')

aln = ProteinAlignment(proteins['BID-201'], proteins['BID-202'])
print(repr(aln))
print(aln.full)

#%%
for protein in proteins[-1:]:
    print(protein.orf)
    for aa in protein.residues:
        if all(aa.codon):
            print(f'\t{aa} <- {aa.codon_str} <- {aa.exons}')

# %%
