#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from sqlalchemy import select

#%%
statement = select(Protein, Transcript).join(Protein.orf).join(ORF.transcript).filter(Transcript.name.contains('MAPK1'))
result = db_session.execute(statement).all()

proteins = [row[Protein] for row in result]
transcripts = [row[Transcript] for row in result]

#%%
for protein in proteins[-1:]:
    print(protein.orf)
    for aa in protein.residues:
        if all(aa.codon):
            print(f'\t{aa} <- {aa.codon} <- {aa.exons}')
# %%
for transcript in transcripts[-1:]:
    print(transcript)
    # for nt in transcript.nucleotides:
    #     if nt.amino_acid:
    #         print(f'\t{nt} -> {nt.amino_acid}')
    for exon in transcript.exons:
        print(f'\t{exon}')
        for nt in exon.nucleotides:
            assert nt.exon is exon
            if nt.amino_acid:
                print(f'\t\t{nt} -> {nt.amino_acid}')
# %%
