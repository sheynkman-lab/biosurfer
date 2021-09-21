from sqlalchemy import create_engine
from sqlalchemy.orm import scoped_session, sessionmaker

DB_GENCODE = 'sqlite:///gencode.sqlite3'
DB_BONE = 'sqlite:///bone.sqlite3'
DB_MEMORY = 'sqlite://'


default_sessiomaker = scoped_session(sessionmaker(autocommit=False, autoflush=False))


class Database:
    registry = {}

    def __new__(cls, path):
        if path in Database.registry:
            return Database.registry[path]
        else:
            obj = super().__new__(cls)
            Database.registry[path] = obj
            return obj

    def __init__(self, path, sessionmaker=None):
        self.path = path
        self._engine = None
        if sessionmaker is None:
            self.sessionmaker = default_sessiomaker
        else:
            self.sessionmaker = sessionmaker
    
    @property
    def engine(self):
        if self._engine is None:
            self._engine = create_engine(self.path)
        return self._engine

    def get_session(self):
        return self.sessionmaker(bind=self.engine)
