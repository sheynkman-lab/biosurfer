import pytest
from Bio.Seq import Seq
from biosurfer.core.models import Transcript, Exon, ORF, Protein

@pytest.fixture(scope='module')
def dummy_utr5_seq():
    return ('GTTTAAGGGATACTTGGTCCATAATGTGCGATACCGACTCATTCAGCCTATAAAGTCTGT'
            'GCAGACTAGAGGCCATTGGAAGACCGAGCCCAACGTCGTA')

@pytest.fixture(scope='module')
def dummy_orf_seq():
    return ('ATGGGCTTTCGCTGCTATCAGATTATTAGCCGCAGCAACGCGATTCTGGTGGTGTGCGGC'
            'CCGTTTCTGGAAGAAGCGGTGGTGTTTAGCCGCGGCAGCCTGTATGATCGCGCGACCTTT'
            'TTTTATAGCGCGGGCATTGCGAGCTTTGAACTGTTTGAACCGCGCGAACATCCGGAACGC'
            'TTTCTGCTGAGCTGCCGCGTGCATATTGTGATGAGCGTGGCGAAAGGCAAAGGCGCGCGC'
            'GGCCGCGCGCAGGTGCAGGTGCCGCTGAGCCCGCGCCTGGCGAGCGGCATTCGCGGCATG'
            'CTGCCGAACCTGGTGACCGCGTGGGCGCGCCGCCCGTGCGCGGATGAAGAATATCAGGGC'
            'GGCAACTATGCGGGCAGCGAAGGCCTGCGCAGCCTGTTTACCAAATATTGGGGCTGCCGC'
            'CTGTATCTGCCGATTCATCAGCGCCTGGCGAACGATGATCCGACCTATAGCTGCCATCCG'
            'CGCGATCTGCAGCAGTTTATTCCGGGCAGCCTGCGCGGCCAGGGCGAAGCGAGCAGCGCG'
            'GCGAGCGGCCCGCATGATAACCGCTGCGCGCGCACCGGCAAACTGATTGCGTGGCTGTGA')

@pytest.fixture(scope='module')
def dummy_utr3_seq():
    return ('AGATCTTTTCACTCCAAGACGTTCTTGTCTGTCATCTGCGCGATATCATAATTGGAAAGC'
            'TGCCGCGCCGGTTGTTGTTGCGCGAACGGTGCCTCAGCTGCGATCCGGGCTCGCTTGTCG'
            'CGTCCTTCCAAAGCTTGTAGTGGAAGTAGCTCGCCCGAGACAGAAGGGCACCTTATGGGG'
            'TGATCTTTGCCTTAGTGAATAACTGATCCACGCTAACCGACCTGGAGTCTATTGAGTGTC'
            'CGACTAGAATATCTAAGGCTCGCAGGCCGCATCCATTGTAACTACCGGACGCGCCGTCGA')

@pytest.fixture(scope='module')
def dummy_exon_ranges():
    return [
        (1, 151),
        (152, 251),
        (252, 386),
        (387, 493),
        (494, 649),
        (650, 1000)
    ]

@pytest.fixture(scope='module')
def dummy_transcript(dummy_utr5_seq, dummy_orf_seq, dummy_utr3_seq, dummy_exon_ranges):
    transcript = Transcript(
        sequence =  dummy_utr5_seq + dummy_orf_seq + dummy_utr3_seq,
        name = 'DUMMY-1'
    )
    for i, (start, stop) in enumerate(dummy_exon_ranges, start=1):
        Exon(position=i, transcript_start=start, transcript_stop=stop, transcript=transcript)
    return transcript

@pytest.fixture(scope='module')
def dummy_protein(dummy_orf_seq):
    protein_seq = Seq(dummy_orf_seq).translate(to_stop=True)
    return Protein(
        sequence =  str(protein_seq),
    )

@pytest.fixture(scope='module')
def dummy_orf(dummy_transcript, dummy_protein, dummy_orf_seq):
    return ORF(
        transcript = dummy_transcript,
        protein = dummy_protein,
        transcript_start = 1,
        transcript_stop = len(dummy_orf_seq)
    )
