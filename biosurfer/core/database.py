from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base, declared_attr, has_inherited_table
from sqlalchemy.orm import scoped_session, sessionmaker

working_dir = '/home/redox/sheynkman-lab/biosurfer/biosurfer/core'
# db_path = f'sqlite:///{working_dir}/gencode.sqlite3'
db_path = f'sqlite:///{working_dir}/test.sqlite3'
# db_path = f'sqlite:///{working_dir}/wtc11.sqlite3'
# db_path = 'sqlite://'


engine = create_engine(db_path, convert_unicode=True)
db_session = scoped_session(sessionmaker(autocommit=False,
                                         autoflush=False,
                                         bind=engine))

class Base:
    @declared_attr
    def __tablename__(cls):
        if has_inherited_table(cls):
            return None
        return cls.__name__.lower()
Base = declarative_base(cls=Base)
Base.query = db_session.query_property()
