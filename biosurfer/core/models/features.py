import json
from functools import cached_property
from typing import TYPE_CHECKING, List

from biosurfer.core.constants import FeatureType
from biosurfer.core.constants import \
    TranscriptLevelAlignmentCategory as TranscriptAlignCat
from biosurfer.core.helpers import run_length_decode
from biosurfer.core.models.base import AccessionMixin, Base, NameMixin, TablenameMixin
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship
from sqlalchemy.sql.schema import Column, ForeignKey, Table, UniqueConstraint
from sqlalchemy.sql.sqltypes import Boolean, Enum, Integer, PickleType, String

if TYPE_CHECKING:
    from biosurfer.core.alignments import FeatureAlignment
    from biosurfer.core.models.biomolecules import Residue

# feature_base_table = Table(
#     'proteinfeature', Base.metadata,
#     Column('type', Enum(FeatureType)),
#     Column('accession', String, primary_key=True, index=True),
#     Column('name', String),
#     Column('description', String)
# )


# feature_mapping_table = Table(
#     'proteinfeature_mapping', Base.metadata,
#     Column('feature_id', String, ForeignKey('proteinfeature.accession'), primary_key=True),
#     Column('protein_id', String, ForeignKey('protein.accession'), primary_key=True),
#     Column('protein_start', Integer, primary_key=True),
#     Column('protein_stop', Integer, primary_key=True),
# )


class Feature(Base, TablenameMixin, NameMixin, AccessionMixin):
    type = Column(Enum(FeatureType))
    description = Column(String)

    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': FeatureType.NONE
    }


class Domain(Feature):
    __tablename__ = None
    __mapper_args__ = {
        'polymorphic_identity': FeatureType.DOMAIN
    }


class ProteinFeature(Base, TablenameMixin):
    id = Column(Integer, primary_key=True, autoincrement=True)
    feature_id = Column(String, ForeignKey('feature.accession'), nullable=False)
    protein_id = Column(String, ForeignKey('protein.accession'), nullable=False)
    protein_start = Column(Integer, nullable=False)
    protein_stop = Column(Integer, nullable=False)
    reference = Column(Boolean, nullable=False)

    feature = relationship('Feature', uselist=False, lazy='selectin')
    protein = relationship('Protein', back_populates='features', uselist=False)

    __table_args__ = (UniqueConstraint(feature_id, protein_id, protein_start, protein_stop),)

    __mapper_args__ = {
        'polymorphic_on': reference,
        'polymorphic_identity': True
    }
    
    def __repr__(self):
        return f'{self.protein}:{self.name}({self.protein_start}-{self.protein_stop})'

    @property
    def type(self) -> FeatureType:
        return self.feature.type if self.feature else FeatureType.NONE

    @property
    def name(self) -> str:
        return self.feature.name if self.feature else None
    
    @property
    def description(self) -> str:
        return self.feature.description if self.feature else None

    @hybrid_property
    def length(self):
        return self.protein_stop - self.protein_start + 1

    @property
    def sequence(self) -> str:
        return self.protein.sequence[self.protein_start-1:self.protein_stop]
    
    @property
    def residues(self) -> List['Residue']:
        return self.protein.residues[self.protein_start-1:self.protein_stop]


class ProjectedFeature(ProteinFeature, TablenameMixin):
    __tablename__ = None
    anchor_id = Column(Integer, ForeignKey('proteinfeature.id'))
    anchor = relationship('ProteinFeature', foreign_keys=[anchor_id], uselist=False)
    _differences = Column(String)  # run-length encoding of FeatureAlignment with anchor feature

    __mapper_args__ = {
        'polymorphic_identity': False
    }

    @cached_property
    def altered_residues(self):
        alt_res = []
        i = 0
        for token in self._differences.split(','):
            char = token[-1]
            length = int(token[:-1])
            if TranscriptAlignCat(char) not in {TranscriptAlignCat.MATCH, TranscriptAlignCat.EDGE_MATCH}:
                alt_res.extend(self.residues[i:i+length])
            i += length
        return alt_res

    # def __init__(self, feature_alignment: 'FeatureAlignment'):
    #     self.alignment = feature_alignment
    #     self.anchor = feature_alignment.proteinfeature
    #     self.feature = feature_alignment.proteinfeature.feature
    #     self.protein = feature_alignment.other
        
    #     residues = feature_alignment.other_residues
    #     self.protein_start = residues[0].position
    #     self.protein_stop = residues[-1].position

    #     self.altered_residues = [res_aln.other for res_aln in feature_alignment if res_aln.category not in {TranscriptAlignCat.MATCH, TranscriptAlignCat.EDGE_MATCH, TranscriptAlignCat.DELETION}]

    def __repr__(self):
        return super().__repr__() + '*'
