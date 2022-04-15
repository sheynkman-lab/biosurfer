from itertools import product
from more_itertools import only, partition
from statistics import median
from typing import TYPE_CHECKING, Iterable, List
from biosurfer.core.alignments import TranscriptAlignment, CodonAlignment, ProteinAlignment

from biosurfer.core.constants import (AnnotationFlag,
                                      SequenceAlignmentCategory,
                                      ProteinRegion, Strand)
from biosurfer.core.models.nonpersistent import Junction
from biosurfer.core.splice_events import ExonBypassEvent

if TYPE_CHECKING:
    from biosurfer.core.splice_events import BasicTranscriptEvent
    from biosurfer.core.alignments import CodonAlignmentBlock, ProteinAlignmentBlock
    from biosurfer.core.models.nonpersistent import Junction, Protein, Transcript


def split_transcripts_on_junction_usage(junction: 'Junction', transcripts: Iterable['Transcript']):
    return partition(lambda transcript: junction in transcript.junctions, transcripts)

def pairwise_align_protein_isoforms(protein_group_a: Iterable['Protein'], protein_group_b: Iterable['Protein']):
    return [
        (
            TranscriptAlignment.from_transcripts(anchor.transcript, other.transcript),
            CodonAlignment.from_proteins(anchor, other),
            ProteinAlignment.from_proteins(anchor, other)
        ) for anchor, other in product(protein_group_a, protein_group_b)
    ]

def get_transcript_events_associated_with_junction(junction: 'Junction', tx_aln: 'TranscriptAlignment'):
    return [
        event for event in tx_aln.basic_events
        if (
            junction in getattr(event, 'anchor_junctions', ()) + getattr(event, 'other_junctions', ())
            or isinstance(event, ExonBypassEvent) and (event.exon & Junction(junction.donor - 1, junction.acceptor + 1))
        )
    ]

def get_cblocks_attributed_to_transcript_event(tx_event: 'BasicTranscriptEvent', cd_aln: 'CodonAlignment'):
    tx_aln = TranscriptAlignment.from_transcripts(cd_aln.anchor.transcript, cd_aln.other.transcript)
    tblock = tx_aln.event_to_block.get(tx_event)
    return cd_aln.tblock_to_cblocks.get(tblock, ())

# def junction_has_drastic_effect_in_pair(
#     junction: 'Junction' = None,
#     anchor: 'Protein' = None,
#     other: 'Protein' = None,
#     pblocks: Iterable['ProteinAlignmentBlock'] = None,
#     threshold_delta_length: int = None) -> bool:

#     if pblocks is None:
#         pblocks = get_pblocks_attributed_to_junction(junction, [Alignment(anchor, other)])
#     else:
#         anchor = pblocks[0].anchor
#         other = pblocks[0].other
    
#     if len(pblocks) == 0:
#         return False
#     if threshold_delta_length is None:
#         threshold_delta_length = max(anchor.length, other.length) * 2 // 5
#     # FIXME: this will return a false positive if junction-related pblock is deletion, but another pblock is a similar-length insertion
#     return (anchor.orf.nmd ^ other.orf.nmd or
#             abs(median(pblock.delta_length for pblock in pblocks)) >= abs(threshold_delta_length))

# def get_event_counts(pblocks: Iterable['ProteinAlignmentBlock']):
#     events = [event for event in AnnotationFlag.__members__.values() if event is not AnnotationFlag.NONE]
#     counts = {event: sum(event & pblock.flags == event for pblock in pblocks) for event in events}
#     return counts
