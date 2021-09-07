from biosurfer.core.constants import Nucleobase
from biosurfer.core.models import Transcript, scoped_session
from biosurfer.core.variant import Variant,  VariantBuilder, VariantTranscript
import unittest

class TestVariantBuilder(unittest.TestCase):
    def test_variant_transcript(self):
        transcript = Transcript.from_accession('ENST00000275493.7')
        variant = Variant(
            chromosome = transcript.chromosome,
            chromosome_position = 55019018,
            reference_nucleotide_sequence = 'G',
            variant_nucleotide_sequence = 'C' 
        )
        variant_builder = VariantBuilder(transcript, [variant])
        variant_transcript = variant_builder.variant_transcript
        assert variant_transcript.get_nucleotide_from_coordinate(55019018).base == Nucleobase('C') 

        
        

        
