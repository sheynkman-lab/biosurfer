import pytest
from biosurfer.core.database import Database
from biosurfer.core.models.base import Base
from biosurfer.core.helpers import get_ids_from_gencode_fasta
from hypothesis import settings
from pathlib import Path

data_dir = Path(__file__).parent.parent/'data'/'gencode'

def pytest_sessionstart(session):
    settings.register_profile('default', deadline=500)
    settings.load_profile('default')

@pytest.fixture(scope='session')
def database_path(tmp_path_factory):
    return tmp_path_factory.mktemp('database')/('test.sqlite3')

@pytest.fixture(scope='session')
def database(database_path):
    db_url = f'sqlite:///{database_path}'
    db = Database(url=db_url)
    # db.engine.echo = True
    Base.metadata.drop_all(db.engine)
    Base.metadata.create_all(db.engine)
    db.load_gencode_gtf(data_dir/'gencode.v38.toy.gtf', overwrite=True)
    db.load_transcript_fasta(
        data_dir/'gencode.v38.toy.transcripts.fa',
        id_extractor=get_ids_from_gencode_fasta,
    )
    db.load_translation_fasta(
        data_dir/'gencode.v38.toy.translations.fa',
        id_extractor=get_ids_from_gencode_fasta,
    )
    db.load_domains(data_dir/'pfamA.tsv')
    db.load_patterns(data_dir/'prosite.dat')
    db.load_feature_mappings(data_dir/'grch38-protein-features.tsv')
    db.project_feature_mappings()
    yield db
    db.engine.dispose()

@pytest.fixture(scope='session')
def session(database):
    session = database.get_session()
    print(session.get_bind())
    yield session
    session.close()

