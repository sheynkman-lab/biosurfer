from typing import Iterable, Type
from sqlalchemy import Column, String, select
from sqlalchemy.ext.declarative import (DeclarativeMeta, declarative_base,
                                        has_inherited_table)
from sqlalchemy.orm.exc import NoResultFound


class BaseMeta(DeclarativeMeta):
    _session = None

    def __init__(cls, classname, bases, dict_, **kwargs):
        cls.__tablename__ = cls.__name__.lower() if not has_inherited_table(cls) else None
        DeclarativeMeta.__init__(cls, classname, bases, dict_, **kwargs)

    @property
    def session(cls):
        if BaseMeta._session is None:
            raise AttributeError('Base.session has not been set')
        return BaseMeta._session
    
    @session.setter
    def session(cls, session):
        BaseMeta._session = session


Base = declarative_base(metaclass=BaseMeta)


class NameMixin:
    name = Column(String, index=True)

    @classmethod
    def from_name(cls: Type['Base'], name: str, unique: bool = True):
        statement = select(cls).where(cls.name == name)
        result = cls.session.execute(statement).scalars()
        if unique:
            try:
                return result.one()
            except NoResultFound:
                return None
        else:
            return result.all()

    @classmethod
    def from_names(cls: Type['Base'], names: Iterable[str]):
        statement = select(cls).where(cls.name.in_(names))
        return {inst.name: inst for inst in cls.session.execute(statement).scalars()}


class AccessionMixin:
    accession = Column(String, primary_key=True, index=True)

    @classmethod
    def from_accession(cls: Type['Base'], accession: str):
        statement = select(cls).where(cls.accession == accession)
        result = cls.session.execute(statement).scalars()
        try:
            return result.one()
        except NoResultFound:
            return None
    
    @classmethod
    def from_accessions(cls: Type['Base'], accessions: Iterable[str]):
        statement = select(cls).where(cls.accession.in_(accessions))
        return {inst.name: inst for inst in cls.session.execute(statement).scalars()}

