from itertools import chain

def test_exon_lengths_correct(transcript):
    assert sum(exon.length for exon in transcript.exons) == transcript.length

def test_exon_sequences_correct(transcript):
    assert ''.join(exon.sequence for exon in transcript.exons) == transcript.sequence

def test_nucleotide_has_residue(transcript):
    orf = transcript.primary_orf
    assert list(chain.from_iterable((res, res, res) for res in orf.protein.residues)) == [nt.residue for nt in orf.nucleotides]

def test_residue_has_nucleotide(transcript):
    orf = transcript.primary_orf
    assert list(chain.from_iterable(res.codon for res in orf.protein.residues)) == orf.nucleotides
