from biosurfer.core.alignments import TranscriptBasedAlignment as Alignment
from biosurfer.core.models import ProteinFeature, Transcript

def test_alignment_full(session, alignment_case):
    Transcript.session = session
    print(session.get_bind())
    anchor = alignment_case['anchor']
    other = alignment_case['other']
    expected = alignment_case['full']
    txs = Transcript.from_names((anchor, other))
    aln = Alignment(txs[anchor].protein, txs[other].protein)
    assert aln.full == expected

def test_feature_alignment(session, feature_case):
    Transcript.session = session
    print(session.get_bind())
    anchor = feature_case['anchor']
    other = feature_case['other']
    txs = Transcript.from_names((anchor, other))
    aln = Alignment(txs[anchor].protein, txs[other].protein)
    for feat_aln, expected in zip(aln.feature_alignments, feature_case['expected']):
        assert feat_aln.full == expected
