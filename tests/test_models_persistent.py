from itertools import chain
from operator import attrgetter

from biosurfer.core.constants import AminoAcid
from biosurfer.core.models import Gene
from sqlalchemy import select

def pytest_generate_tests(metafunc):
    db_session = metafunc.config.db_session
    example_genes = db_session.execute(select(Gene)).scalars().all()
    example_transcripts = list(chain.from_iterable(gene.transcripts for gene in example_genes))
    example_features = [feature for tx in example_transcripts for feature in tx.protein.features if tx.protein and tx.protein.features]

    if 'gene' in metafunc.fixturenames:
        metafunc.parametrize('gene', example_genes, ids=attrgetter('name'))
    if 'transcript' in metafunc.fixturenames:
        metafunc.parametrize('transcript', example_transcripts, ids=attrgetter('name'))
    if 'feature' in metafunc.fixturenames:
        metafunc.parametrize('feature', example_features, ids=str)

def test_exon_lengths_correct(transcript):
    assert sum(exon.length for exon in transcript.exons) == transcript.length

def test_exon_sequences_correct(transcript):
    assert ''.join(exon.sequence for exon in transcript.exons) == transcript.sequence

def test_nucleotide_has_exon(transcript):
    assert all(nt.exon is exon for exon in transcript.exons for nt in exon.nucleotides)

def test_nucleotide_has_residue(transcript):
    orf = transcript.primary_orf
    assert list(chain.from_iterable((res, res, res) for res in orf.protein.residues)) == [nt.residue for nt in orf.nucleotides]

def test_residue_has_nucleotide(transcript):
    orf = transcript.primary_orf
    assert list(chain.from_iterable(res.codon for res in orf.protein.residues)) == orf.nucleotides

def test_feature_residues_correct(feature):
    assert [(res.position, res.amino_acid) for res in feature.residues] == [(feature.protein_start + i, AminoAcid(aa)) for i, aa in enumerate(feature.sequence)]
