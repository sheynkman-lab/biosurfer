from typing import Iterable, Type
from sqlalchemy import Column, String, select
from sqlalchemy.ext.declarative import declarative_base, declared_attr
from sqlalchemy.orm.exc import NoResultFound

Base = declarative_base()


class TablenameMixin:
    @declared_attr
    def __tablename__(cls: Type['Base']):
        return cls.__name__.lower()


class NameMixin:
    name = Column(String, index=True)

    @classmethod
    def from_name(cls: Type['Base'], session, name: str, unique: bool = True):
        statement = select(cls).where(cls.name == name)
        result = session.execute(statement).scalars()
        if unique:
            try:
                return result.one()
            except NoResultFound:
                return None
        else:
            return result.all()

    @classmethod
    def from_names(cls: Type['Base'], session, names: Iterable[str]):
        statement = select(cls).where(cls.name.in_(names))
        return {inst.name: inst for inst in session.execute(statement).scalars()}


class AccessionMixin:
    accession = Column(String, primary_key=True, index=True)

    @classmethod
    def from_accession(cls: Type['Base'], session, accession: str):
        statement = select(cls).where(cls.accession == accession)
        result = session.execute(statement).scalars()
        try:
            return result.one()
        except NoResultFound:
            return None
    
    @classmethod
    def from_accessions(cls: Type['Base'], session, accessions: Iterable[str]):
        statement = select(cls).where(cls.accession.in_(accessions))
        return {inst.name: inst for inst in session.execute(statement).scalars()}
