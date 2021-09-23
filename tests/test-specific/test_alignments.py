from biosurfer.core.alignments import TranscriptBasedAlignment as Alignment
from biosurfer.core.models import Transcript

def test_alignment_full(session, alignment_case):
    anchor = alignment_case['anchor']
    other = alignment_case['other']
    expected = alignment_case['full']
    Transcript.session = session
    txs = Transcript.from_names((anchor, other))
    aln = Alignment(txs[anchor].protein, txs[other].protein)
    assert aln.full == expected
