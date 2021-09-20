from biosurfer.core.models import Protein, Residue

class ProteinFeature:
    def __init__(self, protein: 'Protein', protein_start: int, protein_stop: int):
        self.protein = protein
        self.protein_start = protein_start
        self.protein_stop = protein_stop
    
    @property
    def sequence(self):
        return self.protein.sequence[self.protein_start-1:self.protein_stop]
    
    @property
    def residues(self):
        return self.protein.residues[self.protein_start-1:self.protein_stop]
    