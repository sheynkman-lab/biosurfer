import pytest

from biosurfer.core.constants import AminoAcid
from biosurfer.core.features import ProteinFeature

@pytest.fixture
def dummy_feature(dummy_protein):
    return ProteinFeature(dummy_protein, protein_start=45, protein_stop=60)

sequence = 'GIASFELFEPREHPER'

def test_feature_relationships(dummy_feature, dummy_protein):
    assert dummy_feature.protein is dummy_protein

def test_feature_sequence(dummy_feature):
    assert dummy_feature.sequence == sequence

def test_feature_residues(dummy_feature):
    assert [(res.position, res.amino_acid) for res in dummy_feature.residues] == [(dummy_feature.protein_start + i, AminoAcid(aa)) for i, aa in enumerate(sequence)]
