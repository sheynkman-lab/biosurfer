import pytest
from biosurfer.core.alignments import Alignment as Alignment
from biosurfer.core.models.biomolecules import Transcript

def test_alignment_full(session, alignment_case):
    print(session.get_bind())
    anchor = alignment_case['anchor']
    other = alignment_case['other']
    expected = alignment_case['full']
    txs = Transcript.from_names(session, (anchor, other))
    aln = Alignment(txs[anchor].protein, txs[other].protein)
    assert aln.full == expected

@pytest.mark.skip(reason='feature cases may be outdated')
def test_feature_alignment(session, feature_case):
    print(session.get_bind())
    anchor = feature_case['anchor']
    other = feature_case['other']
    txs = Transcript.from_names(session, (anchor, other))
    aln = Alignment(txs[anchor].protein, txs[other].protein)
    for feat_aln, expected in zip((aln.project_feature(feat)[1] for feat in aln.anchor.features), feature_case['expected']):
        assert feat_aln.full == expected
