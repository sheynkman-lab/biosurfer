import pytest
from biosurfer.core.database import DB_MEMORY, Database
from biosurfer.core.helpers import get_ids_from_gencode_fasta

def pytest_configure(config):
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
    db.load_domain_mappings(
        '../data/gencode/2019-07-04_HMMER_domain_mappings_to_GS_fasta_file.txt',
        '../data/gencode/pfam_a_names.tsv'
    )
    config.database = db
    config.db_session = db.get_session()

def pytest_unconfigure(config):
    config.db_session.close()
    config.database.engine.dispose()

@pytest.fixture(scope='session')
def database(request):
    return request.config.database

@pytest.fixture(scope='function')
def session(request):
    return request.config.db_session
