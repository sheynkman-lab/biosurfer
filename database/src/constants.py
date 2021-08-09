from enum import Enum


class StringEnum(Enum):
    def __str__(self):
        return self.value


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


class ProteinLevelAlignmentCategory(StringEnum):
    MATCH = 'M'
    INSERTION = 'I'
    DELETION = 'D'
    SUBSTITUTION = 'S'
    OTHER = 'X'


class TranscriptLevelAlignmentCategory(StringEnum):
    MATCH = 'm'
    INSERTION = 'i'
    DELETION = 'd'
    FRAMESHIFT = 'f'
    EDGE_MISMATCH = 'e'
    OTHER = 'x'