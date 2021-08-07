#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from alignments import ProteinAlignment
from sqlalchemy import select
from itertools import combinations

#%%
statement = (select(Protein, Transcript).join(Protein.orf).join(ORF.transcript)
            .filter(Transcript.name.contains('TANGO2')))
result = db_session.execute(statement).all()

proteins = [row[Protein] for row in result]
transcripts = [row[Transcript] for row in result]

#%%
for anchor, other in combinations(proteins, 2):
    aln = ProteinAlignment(anchor, other)
    print(repr(aln))
    for res_aln in aln.chain:
        if res_aln.ttype.value == 'f':
            print(f'\t{res_aln}, {res_aln.anchor.codon_str}|{res_aln.other.codon_str}')

#%%
for protein in proteins[-1:]:
    print(protein.orf)
    for aa in protein.residues:
        if all(aa.codon):
            print(f'\t{aa} <- {aa.codon_str} <- {aa.exons}')
# %%
for transcript in transcripts[-1:]:
    print(transcript)
    for exon in transcript.exons:
        print(f'\t{exon}')
        for nt in exon.nucleotides:
            assert nt.exon is exon
            if nt.residue:
                print(f'\t\t{nt} -> {nt.residue}')
# %%
