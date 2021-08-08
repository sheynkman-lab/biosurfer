#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from alignments import ProteinAlignment
from alignments import TranscriptLevelEvent as TLE
from sqlalchemy import select
from itertools import combinations, islice

def get_gene_protein_isoforms(gene_name):
    gene = db_session.execute(select(Gene).filter(Gene.name == gene_name)).one()[Gene]
    return [transcript.orfs[0].protein for transcript in gene.transcripts]

#%%
proteins = get_gene_protein_isoforms('TANGO2')

for anchor, other in combinations(proteins, 2):
    aln = ProteinAlignment(anchor, other)
    print(repr(aln))
    for res_aln in aln.chain:
        if res_aln.ttype in (TLE.FRAMESHIFT, TLE.SPLIT):
            print(f'\t{res_aln}, {res_aln.anchor.codon_str}|{res_aln.other.codon_str}')

#%%
proteins = get_gene_protein_isoforms('RBFOX2')

for anchor, other in combinations(proteins, 2):
    aln = ProteinAlignment(anchor, other)
    print(repr(aln))
    for res_aln in aln.chain:
        if res_aln.ttype in (TLE.FRAMESHIFT, TLE.SPLIT):
            print(f'\t{res_aln}, {res_aln.anchor.codon_str}|{res_aln.other.codon_str}')

#%%
for protein in proteins[-1:]:
    print(protein.orf)
    for aa in protein.residues:
        if all(aa.codon):
            print(f'\t{aa} <- {aa.codon_str} <- {aa.exons}')

# %%
