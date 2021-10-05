import pytest
from biosurfer.core.database import Database
from biosurfer.core.models import Gene, Transcript

def test_database_registry(database, database_path):
    assert Database(url=f'sqlite:///{database_path}') is database

def test_getting_class_session_before_setting_raises_error():
    with pytest.raises(AttributeError):
        Gene.session

def test_classes_share_session(session):
    Gene.session = session
    assert Gene.session is Transcript.session
