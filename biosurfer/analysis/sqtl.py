from itertools import filterfalse, tee
from typing import Iterable

from biosurfer.core.alignments import pairwise_align_protein_sets
from biosurfer.core.constants import ProteinLevelAlignmentCategory, ProteinRegion, Strand
from biosurfer.core.models import Junction, Transcript


def split_transcripts_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    def contains_both_splice_sites(transcript):
        return (transcript.start <= junction.donor <= transcript.stop and
                transcript.start <= junction.acceptor <= transcript.stop)
    def uses_junction(transcript):
        return junction in transcript.junctions
    tx1, tx2 = tee(filter(contains_both_splice_sites, transcripts))
    transcripts_using = filter(uses_junction, tx1)
    transcripts_not_using = filterfalse(uses_junction, tx2)
    return set(transcripts_using), set(transcripts_not_using)

def get_pblocks_related_to_junction(junction: 'Junction', transcripts: Iterable['Transcript']):
    using, not_using = split_transcripts_on_junction_usage(junction, transcripts)
    pblocks = []
    for aln in pairwise_align_protein_sets((tx.protein for tx in not_using), (tx.protein for tx in using)):
        aln.annotate()
        up_exon, down_exon = aln.other.transcript.get_exons_from_junction(junction)
        def is_related_to_junc(pblock):
            if pblock.region is ProteinRegion.INTERNAL:
                if pblock.category is ProteinLevelAlignmentCategory.DELETION:
                    start = pblock.anchor_residues[0].codon[1].coordinate
                    stop = pblock.anchor_residues[-1].codon[1].coordinate
                else:
                    start = pblock.other_residues[0].codon[1].coordinate
                    stop = pblock.other_residues[-1].codon[1].coordinate
                if junction.strand is Strand.PLUS:
                    return junction.donor <= stop + 1 and start - 1 <= junction.acceptor
                elif junction.strand is Strand.MINUS:
                    return junction.donor >= stop - 1 and start + 1 >= junction.acceptor
            elif pblock.region is ProteinRegion.NTERMINUS:
                return down_exon in pblock.other_exons
            elif pblock.region is ProteinRegion.CTERMINUS:
                return up_exon in pblock.other_exons or down_exon in pblock.other_exons
        pblocks.extend(pblock for pblock in aln.protein_blocks 
                       if pblock.category is not ProteinLevelAlignmentCategory.MATCH and is_related_to_junc(pblock))
    return pblocks, using, not_using
