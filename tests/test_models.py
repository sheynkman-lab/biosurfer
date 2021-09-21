import pytest
from biosurfer.core.models import Transcript, ORF, Protein

@pytest.fixture
def dummy_transcript():
    return Transcript(
        sequence =  'ATGGGCTTTAGATGTTACCAAATCATCTCCCGATCTAATGCAATTCTTGTTGTGTGCGGG'
                    'CCCTTTCTGGAAGAGGCGGTGGTTTTTTCTCGGGGCTCACTGTATGATCGAGCAACCTTT'
                    'TTTTATAGTGCGAAATATGCAGATATTTCCCAGAAGTGCTGGGCGCCCGTCCGGGTGCAT'
                    'TTCCTGCTCTCATGTCGAGTGCATATAGTTATGTCGGTAGCAAAGGGCAAAGGAGCCCGT'
                    'GGCAGAGCTCAAGTACAAGTCCCACTTAGCCCGCGCCTCGCAAGCGGTATTCGGGGAATG'
                    'TTGCCTAATTTGGTCACGGCGTGGGCCAGGCGTCCTTGTGCCGACGAAGAGTATCAGGGG'
                    'GGGAACTACGCCGGGAGTGAAGGGTTACGCTCTTTGTTCACCAAATATTGGGGCTGTCGA'
                    'CTTTATCTACCCATCCACCAGAGACTAGCAAATGACGACCCCACTTATAGTTGCCACCCA'
                    'CGCGACTTGCAACAATTTATTCCTGGTTCACTGCGTGGTCAAGGCGAGGCGTCGTCGGCT'
                    'GCCTCTGGCCCACATGACAATAGATGCGCGCGAACGGGGAAACTAATCGCGTGGCTTTGA',
        name = 'DUMMY_TX'
    )

@pytest.fixture
def dummy_protein():
    return Protein(
        sequence =  'MGFRCYQIISRSNAILVVCGPFLEEAVVFSRGSLYDRATFFYSAKYADISQKCWAPVRVH'
                    'FLLSCRVHIVMSVAKGKGARGRAQVQVPLSPRLASGIRGMLPNLVTAWARRPCADEEYQG'
                    'GNYAGSEGLRSLFTKYWGCRLYLPIHQRLANDDPTYSCHPRDLQQFIPGSLRGQGEASSA'
                    'ASGPHDNRCARTGKLIAWL',
        
    )

@pytest.fixture
def dummy_orf(dummy_transcript, dummy_protein):
    return ORF(
        transcript = dummy_transcript,
        protein = dummy_protein,
        transcript_start = 1,
        transcript_stop = dummy_transcript.length
    )


class TestRelationships:
    def test_transcript_has_orf(dummy_transcript, dummy_orf):
        assert dummy_orf in dummy_transcript.orfs

    def test_orf_has_transcript(dummy_transcript, dummy_orf):
        assert dummy_orf.transcript is dummy_transcript

    def test_orf_has_protein(dummy_orf, dummy_protein):
        assert dummy_orf.protein is dummy_protein

    def test_protein_has_orf(dummy_orf, dummy_protein):
        assert dummy_protein.orf is dummy_orf


