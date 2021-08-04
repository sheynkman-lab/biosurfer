#%%
from database import db_session
from models import Gene, Transcript, Exon, ORF, Protein
from sqlalchemy import select

#%%
statement = select(Protein, Transcript).join(Protein.orf).join(ORF.transcript)
result = db_session.execute(statement).all()

proteins = [row[Protein] for row in result]
transcripts = [row[Transcript] for row in result]

#%%
for protein in proteins[:1]:
    print(protein.orf)
    for aa in protein.amino_acids:
        if all(aa.codon):
            print(f'\t{aa} <- {aa.codon}')
# %%
for transcript in transcripts[:1]:
    print(transcript)
    # for nt in transcript.nucleotides:
    #     if nt.amino_acid:
    #         print(f'\t{nt} -> {nt.amino_acid}')
    for i, exon in enumerate(transcript.exons, start=1):
        print(f'\texon {i}')
        for nt in exon.nucleotides:
            if nt.amino_acid:
                print(f'\t\t{nt} -> {nt.amino_acid}')
# %%
