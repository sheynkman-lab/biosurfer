from pyopenms import ProteaseDigestion, AASequence

from biosurfer.core.models import Protein, Residue

class ProteinDigest:
    def __init__(self, digest = 'Trypsin', missed_cleavages = 0) -> None:
        self.digestion = ProteaseDigestion()
        self.digestion.setEnzyme(digest)
        self.missed_cleavages = missed_cleavages

    def _digest_protein(self, protein):
        protein_sequence = AASequence.fromString(protein.sequence)
        base_peptides = []
        self.digestion.digest(protein_sequence, base_peptides)
        base_peptides_string = [pep.toString() for pep in base_peptides]
        base_peptides = []
        seq_index = 0
        for pep in base_peptides_string:
            residues = []
            for i, _ in enumerate(pep):
                residues.append(protein.residues[seq_index+i])
            peptide = Peptide(protein, residues)
            base_peptides.append(peptide)
            seq_index = seq_index + len(pep)
        return base_peptides
        
    def _extract_missed_cleaves(self, base_peptides):
        peptides_with_misses = []
        protein = base_peptides[0].protein
        for index, peptide in enumerate(base_peptides):
            peptides_with_misses.append(peptide)
            for misses in range(1, self.missed_cleavages + 1):
                if index - misses >= 0:
                    peptides = base_peptides[index - misses:index+1]
                    new_peptide = Peptide(protein, [])
                    for pep in peptides:
                        new_peptide = new_peptide + pep
                    peptides_with_misses.append(new_peptide)
        return peptides_with_misses
    
    def digest_protein(self, protein):
        base_peptides = self._digest_protein(protein)
        peptides_with_misses = self._extract_missed_cleaves(base_peptides)
        return peptides_with_misses


class Peptide:
    def __init__(self, protein, residues) -> None:
        self.protein = protein
        self.residues = residues

    def __add__(self, other):
        if self.protein == other.protein:
            residues = self.residues + other.residues
            peptide = Peptide(self.protein, residues)
            return peptide
        else:
            raise ValueError('Peptides must map to the same protein')
    def __repr__(self) -> str:
        return ''.join([f'{res.amino_acid}' for res in self.residues])

    @property
    def previous_residue(self):
        first_residue = self.residues[0]
        if first_residue.position == 1:
            return '.'
        else :
            return self.protein.residues[first_residue.position - 2]
    @property
    def next_residue(self):
        last_residue = self.residues[-1]
        if last_residue.position == len(self.protein.residues):
            return '.'
        else :
            return self.protein.residues[last_residue.position]