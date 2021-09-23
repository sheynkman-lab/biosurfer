from itertools import chain
from operator import attrgetter

from biosurfer.core.models import Gene, ProteinFeature
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
