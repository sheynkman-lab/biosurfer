from enum import auto

from inscripta.biocantor.location.location_impl import Strand

from helpers import OrderedEnum, StringEnum


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
    GAP = '-'


class UTRType(OrderedEnum):
    FIVE_PRIME = auto()
    THREE_PRIME = auto()

    def __str__(self):
        if self is UTRType.FIVE_PRIME:
            return '5utr'
        elif self is UTRType.THREE_PRIME:
            return '3utr'


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
    FRAMESHIFT = 'f'
    EDGE_MISMATCH = 'e'
    COMPLEX = 'x'
    UNKNOWN = '?'


class APPRIS(OrderedEnum):
    PRINCIPAL = 1
    ALTERNATIVE = 2
    NONE = 3
