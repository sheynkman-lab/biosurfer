from itertools import chain
from operator import attrgetter

from biosurfer.core.models import Gene
from sqlalchemy import select

def pytest_generate_tests(metafunc):
    db_session = metafunc.config.db_session
    example_genes = db_session.execute(select(Gene)).scalars().all()
    example_transcripts = list(chain.from_iterable(gene.transcripts for gene in example_genes))
    if 'gene' in metafunc.fixturenames:
        metafunc.parametrize('gene', example_genes, ids=attrgetter('name'))
    if 'transcript' in metafunc.fixturenames:
        metafunc.parametrize('transcript', example_transcripts, ids=attrgetter('name'))
