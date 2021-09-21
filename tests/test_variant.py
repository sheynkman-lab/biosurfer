from biosurfer.core.constants import Nucleobase
from biosurfer.core.models import Transcript, scoped_session
from biosurfer.core.variant import Variant,  VariantBuilder, VariantTranscript
import unittest



class TestVariantBuilder(unittest.TestCase):
    def setUp(self) -> None:

        return super().setUp()


    def test_single_variant_in_utr_plus(self):
        transcript = Transcript.from_accession('ENST00000275493.7')
        variant = Variant(
            chromosome = transcript.chromosome,
            position = 55019018,
            reference_sequence = 'G',
            variant_sequence = 'A' 
        )

        variant_builder = VariantBuilder(transcript, [variant])
        variant_transcript = VariantTranscript.from_accession(variant_builder.variant_transcript['accession'])
        assert(transcript.sequence != variant_transcript.sequence)
        assert(transcript.orfs[0].protein.sequence == transcript.orfs[0].protein.sequence)
        assert(variant_transcript.get_nucleotide_from_coordinate(55019018).base == Nucleobase('A'))
    
    def test_single_variant_in_coding_region_plus(self):
        transcript = Transcript.from_accession('ENST00000275493.7')
        variant = Variant(
            chromosome = transcript.chromosome,
            position = 55019278,
            reference_sequence = 'A',
            variant_sequence = 'G' 
        )

    def test_single_variant_in_coding_region_negative(self):
        transcript = Transcript.from_accession('ENST00000426406.4')
        variant = Variant(
            chromosome = transcript.chromosome,
            position = 451675,
            reference_sequence = 'G',
            variant_sequence = 'A' 
        )
        variant_builder = VariantBuilder(transcript, [variant])
        variant_transcript = VariantTranscript.from_accession(variant_builder.variant_transcript['accession'])

        assert(transcript.sequence != variant_transcript.sequence)
        assert(transcript.orfs[0].protein.sequence != transcript.orfs[0].protein.sequence)
        assert(variant_transcript.get_nucleotide_from_coordinate(55019018).base == Nucleobase('A'))



        # nuclotide G->A
        # codon GGA ->AGA
        # amino acid : Gly -> Arg



        # variant_transcript = variant_builder.variant_transcript
        # assert variant_transcript.get_nucleotide_from_coordinate(55019018).base == Nucleobase('C') 

        
        

        
#%%
from biosurfer.core.constants import Nucleobase
from biosurfer.core.models import Transcript, scoped_session
from biosurfer.core.variant import Variant,  VariantBuilder, VariantTranscript
transcript = Transcript.from_accession('ENST00000275493.7')
variant = Variant(
    chromosome = transcript.chromosome,
    position = 55019278,
    reference_sequence = 'A',
    variant_sequence = 'G' 
)

variant_builder = VariantBuilder(transcript, [variant])
# print(variant_builder.variant_proteins[0].sequence)
# variant_seq = variant_builder.variant_proteins[0].sequence
# ref_seq = transcript.orfs[0].protein.sequence
# print(transcript.sequence == variant_builder.variant_transcript.sequence)
# variant_protein_sequence = variant_builder.variant_proteins[0].sequence
# ref_protein_sequence = transcript.orfs[0].protein.sequence
# print(variant_protein_sequence == ref_protein_sequence)
# %%

# %%
