from enum import Flag, auto, Enum
from typing import Union

from biosurfer.core.helpers import OrderedEnum, StringEnum


class Strand(OrderedEnum):
    PLUS = auto()
    MINUS = auto()
    UNKNOWN = auto()

    def __str__(self):
        if self is Strand.PLUS:
            return '+'
        elif self is Strand.MINUS:
            return '-'
        else:
            return '?'
    
    @classmethod
    def from_symbol(cls, symbol: str) -> 'Strand':
        if symbol == '+':
            return Strand.PLUS
        elif symbol == '-':
            return Strand.MINUS
        else:
            raise ValueError(f'\'{symbol}\' is not a valid strand')


class Nucleobase(StringEnum):
    ADENINE = 'A'
    CYTOSINE = 'C'
    GUANINE = 'G'
    THYMINE = 'T'
    URACIL = 'U'
    GAP = '-'


class AminoAcid(StringEnum):
    ALANINE = 'A'
    ALA = 'A'    
    ISOLEUCINE = 'I'
    ILE = 'I'    
    LEUCINE = 'L'
    LEU = 'L'
    METHIONINE = 'M'
    MET = 'M'
    VALINE = 'V'
    VAL = 'V'
    PHENYLALANINE = 'F'
    PHE = 'F'
    TRYPTOPHAN = 'W'
    TRP = 'W'
    TYROSINE = 'Y'
    TYR = 'Y'
    ASPARAGINE = 'N'
    ASN = 'N'
    CYSTEINE = 'C'
    CYS = 'C'
    GLUTAMINE = 'Q'
    GLN = 'Q'
    SERINE = 'S'
    SER = 'S'
    THREONINE = 'T'
    THR = 'T'
    ASPARTATE = 'D'
    ASP = 'D'
    GLUTAMATE = 'E'
    GLU = 'E'
    ARGININE = 'R'
    ARG = 'R'
    HISTIDINE = 'H'
    HIS = 'H'
    LYSINE = 'K'
    LYS = 'K'
    GLYCINE = 'G'
    GLY = 'G'
    PROLINE = 'P'
    PRO = 'P'
    SELENOCYSTEINE = 'U'
    SEC = 'U'    
    STOP = '*'  # included for ease of use
    UNKNOWN = 'X'
    GAP = '-'


class ProteinLevelAlignmentCategory(StringEnum):
    MATCH = 'M'
    INSERTION = 'I'
    DELETION = 'D'
    SUBSTITUTION = 'S'
    UNKNOWN = '?'


class TranscriptLevelAlignmentCategory(StringEnum):
    MATCH = 'm'
    INSERTION = 'i'
    DELETION = 'd'
    TRANSLATED = 't'
    UNTRANSLATED = 'u'
    FRAME_AHEAD = 'a'
    FRAME_BEHIND = 'b'
    EDGE = 'e'
    EDGE_MATCH = 'e'
    EDGE_MISMATCH = 'g'
    COMPLEX = 'x'
    UNKNOWN = '?'


AlignmentCategory = Union[ProteinLevelAlignmentCategory, TranscriptLevelAlignmentCategory]


class AnnotationFlag(Flag):
    NONE = 0
    SE = auto()
    IE = auto()
    A5SS = auto()
    A3SS = auto()
    IR = auto()
    IX = auto()
    SIF = auto()

    MXIC = auto()
    UIC = auto()
    DIC = auto()
    TIS = auto()
    UP_TIS = UIC | TIS
    DN_TIS = DIC | TIS
    TSS = auto()

    ACTE = auto()
    UTC = auto()
    DTC = auto()
    EXITC = auto()

    def __str__(self):
        raw = super().__str__()
        return raw.split('.')[1]


class ProteinRegion(OrderedEnum):
    NTERMINUS = auto()
    INTERNAL = auto()
    CTERMINUS = auto()

    def __str__(self):
        if self is ProteinRegion.NTERMINUS:
            return 'Nterm'
        elif self is ProteinRegion.INTERNAL:
            return 'internal'
        elif self is ProteinRegion.CTERMINUS:
            return 'Cterm'


class FeatureType(Enum):
    DOMAIN = auto()
    IDR = auto()
    COILED = auto()
    LCR = auto()
    PTM = auto()
    SIGNALP = auto()
    TRANSMEMBRANE = auto()
    NONE = auto()
    # TODO: add more types


class APPRIS(OrderedEnum):
    NONE = auto()
    ALTERNATIVE = auto()
    PRINCIPAL = auto()


class SQANTI(OrderedEnum):
    FSM = auto()
    ISM = auto()
    NIC = auto()
    NNC = auto()
    OTHER = auto()

    def __str__(self):
        return self.name


START_CODON = 'ATG'
STOP_CODONS = {'TGA', 'TAA', 'TAG'}
