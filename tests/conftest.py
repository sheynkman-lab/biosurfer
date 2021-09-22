from itertools import chain
from operator import attrgetter
import pytest
from biosurfer.core.database import Database, DB_MEMORY
from biosurfer.core.helpers import get_ids_from_gencode_fasta
from biosurfer.core.models import Gene
from sqlalchemy import select

db = Database(DB_MEMORY)
db.load_gencode_gtf('../data/gencode/gencode.v38.toy.gtf', overwrite=True)
db.load_transcript_fasta(
    '../data/gencode/gencode.v38.toy.transcripts.fa',
    id_extractor=get_ids_from_gencode_fasta,
)
db.load_translation_fasta(
    '../data/gencode/gencode.v38.toy.translations.fa',
    id_extractor=get_ids_from_gencode_fasta,
)
db_session = db.get_session()

example_genes = db_session.execute(select(Gene)).scalars().all()
example_transcripts = list(chain.from_iterable(gene.transcripts for gene in example_genes))

@pytest.fixture(scope='session')
def database():
    return db

@pytest.fixture(scope='session')
def session():
    yield db_session
    db_session.close()

def pytest_generate_tests(metafunc):
    if 'gene' in metafunc.fixturenames:
        metafunc.parametrize('gene', example_genes, ids=attrgetter('name'))
    if 'transcript' in metafunc.fixturenames:
        metafunc.parametrize('transcript', example_transcripts, ids=attrgetter('name'))
