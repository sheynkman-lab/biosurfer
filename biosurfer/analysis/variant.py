#%%
from sqlalchemy.sql.expression import desc
from biosurfer.core.models import Protein, Transcript, VariantTranscript
from biosurfer.core.variant import (
    get_transcript_accessions_with_variants,
    get_possible_variants,
    VariantBuilder
)

transcripts_with_possible_variants = get_transcript_accessions_with_variants()
print(len(transcripts_with_possible_variants))
i = 0
for accession in transcripts_with_possible_variants:
    i = i+1
    accession = accession[0]
    if i % 100 ==0:
        print(i)
    transcript = Transcript.from_accession(accession)
    possible_variants = get_possible_variants(transcript)
    variant_builder = VariantBuilder(transcript,possible_variants)


#%%
from biosurfer.core.models import db_session
variant_transcripts = db_session.query(VariantTranscript).all()
# %%
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
protein_sequences = []
for transcript in variant_transcripts:
    protein = transcript.orfs[0].protein
    if protein.sequence is not None:
        pseq = SeqRecord(Seq(protein.sequence), protein.accession, description="",name="")
        protein_sequences.append(pseq)

    
# %%
from Bio import SeqIO
SeqIO.write(protein_sequences, '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/biosurfer/biosurfer/analysis/variant_sequences.fasta','fasta')
# %%
from biosurfer.core.peptides import ProteinDigest
# %%
digest = ProteinDigest(missed_cleavages=0)
peptide_sequences = []
for transcript in variant_transcripts:
    protein = transcript.orfs[0].protein
    peptides = digest.digest_protein(protein)
    peptides = set([str(pep) for pep in peptides])
    ref_transcript = Transcript.from_accession(transcript.reference_transcript_id)
    if len(ref_transcript.orfs) > 0:
        ref_protein = Transcript.from_accession(transcript.reference_transcript_id).orfs[0].protein
        if ref_protein is not None:
            ref_peptides = digest.digest_protein(ref_protein)
            ref_peptides = set([str(pep) for pep in ref_peptides])
            variant_peptides = peptides - ref_peptides
            for var_pep in variant_peptides:
                record = SeqRecord(
                    Seq(var_pep),
                    f'{protein.accession}.{var_pep}',
                    name="",
                    description=""
                )
                peptide_sequences.append(record)
SeqIO.write(peptide_sequences, '/Users/bj8th/Documents/Sheynkman-Lab/GitHub/biosurfer/biosurfer/analysis/variant_peptide_sequences.fasta','fasta')




# %%
