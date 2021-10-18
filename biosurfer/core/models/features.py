from typing import TYPE_CHECKING, List

from biosurfer.core.constants import FeatureType
from biosurfer.core.constants import \
    TranscriptLevelAlignmentCategory as TranscriptAlignCat
from biosurfer.core.models.base import AccessionMixin, Base, NameMixin
from sqlalchemy.ext.hybrid import hybrid_property
from sqlalchemy.orm import relationship
from sqlalchemy.sql.schema import Column, ForeignKey, Table, UniqueConstraint
from sqlalchemy.sql.sqltypes import Boolean, Enum, Integer, String

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


class Feature(Base, NameMixin, AccessionMixin):
    type = Column(Enum(FeatureType))
    description = Column(String)

    __mapper_args__ = {
        'polymorphic_on': type,
        'polymorphic_identity': FeatureType.NONE
    }


class Domain(Feature):
    __mapper_args__ = {
        'polymorphic_identity': FeatureType.DOMAIN
    }


class ProteinFeature(Base):
    __tablename__ = 'proteinfeature'
    id = Column(Integer, primary_key=True, autoincrement=True)
    feature_id = Column(String, ForeignKey('feature.accession'))
    protein_id = Column(String, ForeignKey('protein.accession'))
    protein_start = Column(Integer)
    protein_stop = Column(Integer)
    reference = Column(Boolean, nullable=False)

    feature = relationship('Feature', uselist=False)
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
        return self.feature.type

    @property
    def name(self) -> str:
        return self.feature.name
    
    @property
    def description(self) -> str:
        return self.feature.description

    @hybrid_property
    def length(self):
        return self.protein_stop - self.protein_start + 1

    @property
    def sequence(self) -> str:
        return self.protein.sequence[self.protein_start-1:self.protein_stop]
    
    @property
    def residues(self) -> List['Residue']:
        return self.protein.residues[self.protein_start-1:self.protein_stop]


class ProjectedFeature(ProteinFeature):
    anchor_id = Column(Integer, ForeignKey('proteinfeature.id'))
    anchor = relationship('ProteinFeature', uselist=False)

    __mapper_args__ = {
        'polymorphic_identity': False
    }

    def __init__(self, feature_alignment: 'FeatureAlignment'):
        self.alignment = feature_alignment
        self.anchor = feature_alignment.proteinfeature
        self.feature_id = feature_alignment.proteinfeature.feature_id
        self.protein_id = feature_alignment.other.accession
        
        residues = feature_alignment.other_residues
        self.protein_start = residues[0].position
        self.protein_stop = residues[-1].position

        # TODO: persist this in database
        self.altered_residues = [res_aln.other for res_aln in feature_alignment if res_aln.category not in {TranscriptAlignCat.MATCH, TranscriptAlignCat.EDGE_MATCH, TranscriptAlignCat.DELETION}]

    def __repr__(self):
        return f'{self.anchor.protein}>>' + ProteinFeature.__repr__(self)
