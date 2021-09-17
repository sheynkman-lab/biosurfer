from itertools import filterfalse, tee
from statistics import median
from typing import TYPE_CHECKING, Iterable

from biosurfer.core.alignments import pairwise_align_protein_sets
from biosurfer.core.constants import (AnnotationFlag,
                                      ProteinLevelAlignmentCategory,
                                      ProteinRegion, Strand)
from biosurfer.core.models import Junction, Transcript

if TYPE_CHECKING:
    from biosurfer.core.alignments import (ProteinAlignmentBlock,
                                           TranscriptBasedAlignment)
    from biosurfer.core.models import Protein


def split_transcripts_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    # def contains_both_splice_sites(transcript):
    #     return (transcript.start <= junction.donor <= transcript.stop and
    #             transcript.start <= junction.acceptor <= transcript.stop)
    def uses_junction(transcript):
        return junction in transcript.junctions
    tx1, tx2 = tee(transcripts)
    transcripts_using = filter(uses_junction, tx1)
    transcripts_not_using = filterfalse(uses_junction, tx2)
    return set(transcripts_using), set(transcripts_not_using)

def pairwise_align_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    using, not_using = split_transcripts_on_junction_usage(junction, transcripts)
    alns = pairwise_align_protein_sets((tx.protein for tx in not_using), (tx.protein for tx in using))
    return alns, using, not_using

def get_pblocks_related_to_junction(junction: 'Junction', alns: Iterable['TranscriptBasedAlignment']):
    pblocks = []
    for aln in alns:
        up_exon, down_exon = aln.other.transcript.get_exons_from_junction(junction)
        def is_related_to_junc(pblock):
            result = False
            if pblock.category is ProteinLevelAlignmentCategory.DELETION:
                start = pblock.anchor_residues[0].codon[1].coordinate
                stop = pblock.anchor_residues[-1].codon[1].coordinate
            else:
                start = pblock.other_residues[0].codon[1].coordinate
                stop = pblock.other_residues[-1].codon[1].coordinate
            if junction.strand is Strand.PLUS:
                result = junction.donor <= stop + 1 and start - 1 <= junction.acceptor
            elif junction.strand is Strand.MINUS:
                result = junction.donor >= stop - 1 and start + 1 >= junction.acceptor
            if pblock.region is ProteinRegion.NTERMINUS:
                result |= down_exon in pblock.other_exons
            return result
        pblocks.extend(pblock for pblock in aln.protein_blocks 
                       if pblock.category is not ProteinLevelAlignmentCategory.MATCH and is_related_to_junc(pblock))
    return pblocks

def junction_causes_knockdown_in_pair(
    junction: 'Junction' = None,
    anchor: 'Protein' = None,
    other: 'Protein' = None,
    pblocks: Iterable['ProteinAlignmentBlock'] = None,
    threshold_delta_length: int = None) -> bool:

    if pblocks is None:
        pblocks = get_pblocks_related_to_junction(junction, [TranscriptBasedAlignment(anchor, other)])
    else:
        anchor = pblocks[0].anchor
        other = pblocks[0].other
    
    if len(pblocks) == 0:
        return False
    if threshold_delta_length is None:
        threshold_delta_length = max(anchor.length, other.length) * 2 // 5
    # FIXME: this will return a false positive if junction-related pblock is deletion, but another pblock is a similar-length insertion
    return (anchor.orf.nmd ^ other.orf.nmd or
            abs(median(pblock.delta_length for pblock in pblocks)) >= abs(threshold_delta_length))

def get_event_counts(pblocks: Iterable['ProteinAlignmentBlock']):
    events = [event for event in AnnotationFlag.__members__.values() if event is not AnnotationFlag.NONE]
    counts = {event: sum(event & pblock.flags == event for pblock in pblocks) for event in events}
    return counts
