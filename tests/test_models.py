def test_transcript_has_orf(dummy_transcript, dummy_orf):
    assert dummy_orf in dummy_transcript.orfs

def test_orf_has_transcript(dummy_transcript, dummy_orf):
    assert dummy_orf.transcript is dummy_transcript

def test_orf_has_protein(dummy_orf, dummy_protein):
    assert dummy_orf.protein is dummy_protein

def test_protein_has_orf(dummy_orf, dummy_protein):
    assert dummy_protein.orf is dummy_orf

def test_exon_sequences(dummy_transcript, dummy_utr5_seq, dummy_orf_seq, dummy_utr3_seq, dummy_exon_ranges):
    dummy_transcript_seq = dummy_utr5_seq + dummy_orf_seq + dummy_utr3_seq
    exon_sequences = [dummy_transcript_seq[start-1:stop] for start, stop in dummy_exon_ranges]
    assert [exon.sequence for exon in dummy_transcript.exons] == exon_sequences
