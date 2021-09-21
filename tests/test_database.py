import pytest
from biosurfer.core.database import Database, DB_MEMORY
from biosurfer.core.models import Gene, Transcript

@pytest.fixture
def dummy_database():
    return Database(DB_MEMORY)

def test_database_registry(dummy_database):
    assert Database(DB_MEMORY) is dummy_database

def test_getting_class_session_before_setting_raises_error():
    with pytest.raises(AttributeError):
        Gene.session

def test_classes_share_session(dummy_database):
    Gene.session = dummy_database.get_session()
    assert Gene.session is Transcript.session
