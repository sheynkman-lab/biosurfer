import csv
from abc import abstractmethod
from collections import deque
from collections.abc import Sequence
from functools import cached_property, lru_cache
from itertools import chain, groupby, product
from operator import attrgetter
from typing import Dict, Iterable, List, MutableSequence, Optional, Tuple, Union

from Bio.Align import substitution_matrices
from Bio.pairwise2 import align
from biosurfer.core.constants import AminoAcid, AnnotationFlag
from biosurfer.core.constants import \
    SequenceAlignmentCategory as ProteinAlignCat
from biosurfer.core.constants import ProteinRegion, Strand
from biosurfer.core.constants import \
    CodonAlignmentCategory as TranscriptAlignCat
from biosurfer.core.helpers import run_length_encode
from biosurfer.core.models.biomolecules import Exon, Protein, Residue
from biosurfer.core.models.features import ProjectedFeature, ProteinFeature

PBLOCK_FIELDS = ('anchor', 'other', 'region', 'category', 'delta_length', 'event', 'flags', 'annotation')

MATCH_OR_FRAME = {TranscriptAlignCat.MATCH, TranscriptAlignCat.FRAME_AHEAD, TranscriptAlignCat.FRAME_BEHIND}
FRAMESHIFT = {TranscriptAlignCat.FRAME_AHEAD, TranscriptAlignCat.FRAME_BEHIND}

SUBST_MATRIX = substitution_matrices.load('BLOSUM80')

class GapResidue(Residue):
    def __init__(self, protein: 'Protein', position: int, upstream_exon: Optional['Exon'], downstream_exon: Optional['Exon']):
        super().__init__(protein, position)
        self._aa = AminoAcid.GAP
        self.upstream_exon = upstream_exon
        self.downstream_exon = downstream_exon

    @Residue.codon_str.getter
    def codon_str(self):
        return '-'

    @Residue.exons.getter
    def exons(self):
        return list(filter(None, [self.upstream_exon, self.downstream_exon]))
    
    @Residue.primary_exon.getter
    def primary_exon(self):
        return None


# TODO: need to clean up Alignment class hierarchy
class ResidueAlignment:
    def __init__(self, anchor: 'Residue', other: 'Residue', category: 'TranscriptAlignCat'):
        self.anchor = anchor
        self.other = other
        self.category = category
    
    def __repr__(self):
        return f'{self.anchor}|{self.other}'


class ResidueAlignmentSequence(Sequence[ResidueAlignment]):  # this might break compatibility with 3.8 and earlier?
    @property
    def full(self):
        anchor_str = ''.join(str(res.anchor.amino_acid) for res in self)
        anchor_exon_str = ''.join(str(res.anchor.primary_exon.position % 10) if res.anchor.primary_exon else '-' for res in self)
        other_str = ''.join(str(res.other.amino_acid) for res in self)
        other_exon_str = ''.join(str(res.other.primary_exon.position % 10) if res.other.primary_exon else '-' for res in self)
        cat_str = ''.join(str(res.category) for res in self)
        return '\n'.join([anchor_exon_str, anchor_str, cat_str, other_str, other_exon_str])


class AlignmentBlock(ResidueAlignmentSequence):
    def __init__(self, parent: 'Alignment', start: int, end: int):
        self.parent = parent
        self.start = start  # 0-based, inclusive
        self.end = end  # 0-based, exclusive
        self.length = end - start

    @abstractmethod
    def __repr__(self):
        pass

    def __len__(self):
        return self.length

    def __getitem__(self, index) -> Union[ResidueAlignment, List[ResidueAlignment]]:
        if isinstance(index, slice):
            if index.step:
                raise NotImplementedError('AlignmentBlock does not support slices with steps')
            return [self[i] for i in range(*index.indices(self.length))]
        elif isinstance(index, int):
            return self.parent[self._parent_based_index(index)]
        else:
            raise TypeError(f'{index} is not an int or slice')
    
    def _parent_based_index(self, index):
        if index < -self.length or index >= self.length:
            raise IndexError(f'{index} out of range for {self}')
        return (index % self.length) + self.start
    
    @property
    def anchor(self) -> 'Protein':
        return self[0].anchor.protein

    @property
    def other(self) -> 'Protein':
        return self[0].other.protein

    @cached_property
    def anchor_residues(self) -> List['Residue']:
        return [res_aln.anchor for res_aln in self if not res_aln.anchor.is_gap]

    @cached_property
    def anchor_exons(self):
        exons = {anchor_res.primary_exon for anchor_res in self.anchor_residues}
        exons = exons | {exon for res_aln in self for exon in res_aln.anchor.exons if res_aln.anchor.is_gap}
        return exons

    @cached_property
    def anchor_junctions(self):
        return {res_aln.anchor.junction for res_aln in self if res_aln.anchor.junction}

    @cached_property
    def anchor_sequence(self) -> str:
        return ''.join(str(res.amino_acid) for res in self.anchor_residues)

    @cached_property
    def other_residues(self) -> List['Residue']:
        return [res_aln.other for res_aln in self if not res_aln.other.is_gap]

    @cached_property
    def other_exons(self):
        exons = {other_res.primary_exon for other_res in self.other_residues}
        exons = exons | {exon for res_aln in self for exon in res_aln.other.exons if res_aln.other.is_gap}
        return exons

    @cached_property
    def other_junctions(self):
        return {res_aln.other.junction for res_aln in self if res_aln.other.junction}

    @cached_property
    def other_sequence(self) -> str:
        return ''.join(str(res.amino_acid) for res in self.other_residues)
    
    @property
    def delta_length(self) -> int:
        return len(self.other_residues) - len(self.anchor_residues)
    
    @property
    def similarity(self) -> int:
        return align.globalds(
            self.anchor_sequence if self.anchor_sequence else 'X',
            self.other_sequence if self.other_sequence else 'X',
            SUBST_MATRIX, -11, -1, score_only=True)


class TranscriptAlignmentBlock(AlignmentBlock):
    def __init__(self, parent, position, start, end, category: 'TranscriptAlignCat'):
        super().__init__(parent, start, end)
        self.position = position
        self.category = category
        # These attributes are useful in the annotation code
        self._prev_match_or_frame_tblock = None
        self._next_match_or_frame_tblock = None
        self._annotations = []
        self.flags = AnnotationFlag.NONE
    
    def __repr__(self):
        return f'{self.parent}:tblock{self.position}-{self.category}'


class ProteinAlignmentBlock(AlignmentBlock):
    def __init__(self, parent, position, tblocks: Sequence['TranscriptAlignmentBlock'], category: 'ProteinAlignCat'):
        start = min(tblock.start for tblock in tblocks)
        end = max(tblock.end for tblock in tblocks)
        super().__init__(parent, start, end)
        self.position = position
        self.category = category
        self.transcript_blocks = list(tblocks)
        self._annotations = []
        self.region = ProteinRegion.INTERNAL
        self._event = None

    def __repr__(self):
        return f'{self.parent}:pblock{self.position}-{self.category}'
    
    @property
    def annotation(self):
        out = ', \n'.join(filter(None, [annotation for tblock in self.transcript_blocks for annotation in tblock._annotations] + self._annotations))
        return out if out else None
    
    @property
    def event(self):
        return self._event

    # this is pretty kludgy, will replace when Annotation classes are implemented
    @event.setter
    def event(self, event):
        if self.event is None or 'SIF' in event:
            self._event = event
        elif event in {'FS', 'NMD'}:
            self._event = f'{self._event}-{event}'
        else:
            self._event = 'complex'

    @property
    def flags(self):
        result = AnnotationFlag.NONE
        for tblock in self.transcript_blocks:
            result |= tblock.flags
        return result

    def to_dict(self):
        return {field: getattr(self, field) for field in PBLOCK_FIELDS}


class Alignment(ResidueAlignmentSequence):
    # @lru_cache()
    # def __new__(cls, *args, **kwargs):
    #     return object.__new__(cls)

    def __init__(self, anchor: 'Protein', other: 'Protein'):
        # if anchor.orf.gene is not other.orf.gene:
        #     raise ValueError(f'{anchor} and {other} belong to different genes')
        strand = anchor.orf.transcript.strand
        # if strand is not other.orf.transcript.strand:
        #     raise ValueError(f'{anchor.orf.transcript} and {other.orf.transcript} are on different strands')
        self.anchor = anchor
        self.other = other
        self._chain = rough_alignment(anchor, other, strand)
        refine_alignment(self._chain)
        self.transcript_blocks = get_transcript_blocks(self)
        self.protein_blocks = get_protein_blocks(self)
        self._annotate()
        # self.projected_features = list(filter(None, (self.project_feature(feature) for feature in self.anchor.features)))

    def __repr__(self):
        return f'{self.anchor}|{self.other}'

    def __getitem__(self, index):
        return self._chain[index]
    
    def __len__(self):
        return len(self._chain)
    
    # TODO: use Annotation classes in the future
    # TODO: detect NAGNAG splicing
    def _annotate(self) -> None:
        # DELETE_INSERT = {TranscriptAlignCat.DELETION, TranscriptAlignCat.INSERTION}

        strand = self.anchor.orf.transcript.strand
        anchor_transcript = self.anchor.orf.transcript
        other_transcript = self.other.orf.transcript
        nterminal_pblock = None
        cterminal_pblock = None
        upstream_cterm_tblock = None
        upstream_cterm_res_aln = None
        
        # classify internal splicing events
        for pblock in self.protein_blocks:
            if not nterminal_pblock:
                for res_aln in pblock:
                    if res_aln.anchor.position == 1 or res_aln.other.position == 1:
                        nterminal_pblock = pblock
                        pblock.region = ProteinRegion.NTERMINUS
                        break
            if not cterminal_pblock:
                for res_aln in pblock:
                    if res_aln.anchor.amino_acid is AminoAcid.STOP or res_aln.other.amino_acid is AminoAcid.STOP:
                        cterminal_pblock = pblock
                        pblock.region = ProteinRegion.CTERMINUS
                        for tblock in pblock.transcript_blocks:
                            if res_aln in tblock:
                                upstream_cterm_tblock = tblock
                                break
                        upstream_cterm_res_aln = res_aln
                        break
            
            if pblock.category is ProteinAlignCat.MATCH:
                continue
            
            for tblock in pblock.transcript_blocks:
                if tblock.category is TranscriptAlignCat.DELETION:
                    first_exon = tblock[0].anchor.codon[2].exon
                    last_exon = tblock[-1].anchor.codon[0].exon

                    if tblock._prev_match_or_frame_tblock and tblock._next_match_or_frame_tblock:
                        prev_anchor_exon = tblock._prev_match_or_frame_tblock[-1].anchor.codon[0].exon
                        next_anchor_exon = tblock._next_match_or_frame_tblock[0].anchor.codon[2].exon
                        if prev_anchor_exon is next_anchor_exon:
                            tblock._annotations.append(f'portion of {prev_anchor_exon} intronized')
                            tblock.flags |= AnnotationFlag.IX
                            pblock.event = 'IX'
                        else:
                            e_first = first_exon.position
                            e_last = last_exon.position
                            if prev_anchor_exon is first_exon:
                                tblock._annotations.append(f'{first_exon} shortened by alternative splice donor')
                                tblock.flags |= AnnotationFlag.A5SS
                                pblock.event = 'A5SS-del'
                                e_first += 1
                            if next_anchor_exon is last_exon:
                                tblock._annotations.append(f'{last_exon} shortened by alternative splice acceptor')
                                tblock.flags |= AnnotationFlag.A3SS
                                pblock.event = 'A3SS-del'
                                e_last -= 1
                            first_skipped_exon = anchor_transcript.exons[e_first-1]
                            last_skipped_exon = anchor_transcript.exons[e_last-1]
                            if first_skipped_exon is last_skipped_exon:
                                tblock._annotations.append(f'{first_skipped_exon} skipped')
                                tblock.flags |= AnnotationFlag.SE
                                pblock.event = 'SE'
                            elif e_first < e_last:
                                tblock._annotations.append(f'exons {first_skipped_exon} to {last_skipped_exon} skipped')
                                tblock.flags |= AnnotationFlag.SE
                                pblock.event = 'SE'
                
                elif tblock.category is TranscriptAlignCat.INSERTION:
                    first_exon = tblock[0].other.codon[2].exon
                    last_exon = tblock[-1].other.codon[0].exon

                    if tblock._prev_match_or_frame_tblock and tblock._next_match_or_frame_tblock:
                        prev_anchor_exon = tblock._prev_match_or_frame_tblock[-1].anchor.codon[0].exon
                        next_anchor_exon = tblock._next_match_or_frame_tblock[0].anchor.codon[2].exon
                        prev_other_exon = tblock._prev_match_or_frame_tblock[-1].other.codon[0].exon
                        next_other_exon = tblock._next_match_or_frame_tblock[0].other.codon[2].exon
                        if prev_other_exon is next_other_exon:
                            tblock._annotations.append(f'retained intron between {prev_anchor_exon} and {next_anchor_exon}')
                            tblock.flags |= AnnotationFlag.IR
                            pblock.event = 'IR'
                        else:
                            number_of_included_exons = last_exon.position - first_exon.position + 1
                            if prev_other_exon is first_exon:
                                tblock._annotations.append(f'{prev_anchor_exon} lengthened by alternative splice donor')
                                tblock.flags |= AnnotationFlag.A5SS
                                pblock.event = 'A5SS-ins'
                                number_of_included_exons -= 1
                            if next_other_exon is last_exon:
                                tblock._annotations.append(f'{next_anchor_exon} lengthened by alternative splice acceptor')
                                tblock.flags |= AnnotationFlag.A3SS
                                pblock.event = 'A3SS-ins'
                                number_of_included_exons -= 1
                            if number_of_included_exons == 1:
                                tblock._annotations.append(f'exon included between {prev_anchor_exon} and {next_anchor_exon}')
                                tblock.flags |= AnnotationFlag.IE
                                pblock.event = 'IE'
                            elif number_of_included_exons > 1:
                                tblock._annotations.append(f'{number_of_included_exons} exons included between {prev_anchor_exon} and {next_anchor_exon}')
                                tblock.flags |= AnnotationFlag.IE
                                pblock.event = 'IE'
                
                elif tblock.category in FRAMESHIFT:
                    first_exon = tblock[0].anchor.codon[2].exon
                    last_exon = tblock[-1].anchor.codon[0].exon
                    if first_exon is last_exon:
                        pblock._annotations.append(f'{first_exon} translated in different frame')
                    else:
                        pblock._annotations.append(f'{first_exon} to {last_exon} translated in different frame')
                    tblock.flags |= AnnotationFlag.SIF
                    pblock.event = 'FS'
        
        # classify N-terminal changes (if any)
        if nterminal_pblock and nterminal_pblock.category is not ProteinAlignCat.MATCH:
            upstream_start_codon = tuple(nt.coordinate for nt in self.anchor.residues[0].codon)
            downstream_start_codon = tuple(nt.coordinate for nt in self.other.residues[0].codon)
            upstream_start_transcript, downstream_start_transcript = anchor_transcript, other_transcript
            if strand is Strand.PLUS:
                alternative_tss = anchor_transcript.start != other_transcript.start
                if upstream_start_codon[1] > downstream_start_codon[1]:
                    upstream_start_codon, downstream_start_codon = downstream_start_codon, upstream_start_codon 
                    upstream_start_transcript, downstream_start_transcript = downstream_start_transcript, upstream_start_transcript
            elif strand is Strand.MINUS:
                alternative_tss = anchor_transcript.stop != other_transcript.stop
                if upstream_start_codon[1] < downstream_start_codon[1]:
                    upstream_start_codon, downstream_start_codon = downstream_start_codon, upstream_start_codon 
                    upstream_start_transcript, downstream_start_transcript = downstream_start_transcript, upstream_start_transcript
            
            upstream_start_codon_shared_nts = sum(downstream_start_transcript.contains_coordinate(coord) for coord in upstream_start_codon)
            downstream_start_codon_shared_nts = sum(upstream_start_transcript.contains_coordinate(coord) for coord in downstream_start_codon)
            if upstream_start_codon_shared_nts == 3:
                if downstream_start_codon_shared_nts == 3:
                    # mutually shared start codons
                    if other_transcript is downstream_start_transcript:
                        nterminal_pblock._annotations.append('usage of downstream alternative TIS')  # TODO: indicate anchor exon
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.DN_TIS
                        nterminal_pblock.event = 'dnTIS'
                    else:
                        nterminal_pblock._annotations.append('usage of upstream alternative TIS')  # TODO: indicate anchor exon
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.UP_TIS
                        nterminal_pblock.event = 'upTIS'
                else:
                    # shared upstream start, exclusive downstream start
                    if other_transcript is downstream_start_transcript:
                        nterminal_pblock._annotations.append('usage of downstream alternative TIS revealed by splicing')  # TODO: indicate surrounding anchor exons
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.DN_TIS
                        nterminal_pblock.event = 'dnTIS'
                    else:
                        nterminal_pblock._annotations.append('downstream start codon spliced out leading to usage of upstream start codon')  # TODO: indicate anchor exon
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.UIC
                        nterminal_pblock.event = 'UIC-splice'
            else:
                if downstream_start_codon_shared_nts == 3:
                    # exclusive upstream start, shared downstream start
                    if strand is Strand.PLUS:
                        caused_by_alt_tss = upstream_start_codon[0] < downstream_start_transcript.start
                    elif strand is Strand.MINUS:
                        caused_by_alt_tss = upstream_start_codon[0] > downstream_start_transcript.stop
                    if other_transcript is downstream_start_transcript:
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.DIC
                        if caused_by_alt_tss:
                            for tblock in nterminal_pblock.transcript_blocks:
                                tblock.flags |= AnnotationFlag.TSS
                            nterminal_pblock._annotations.append('alternative TSS leading to usage of downstream start codon')  # TODO: indicate anchor exon
                            nterminal_pblock.event = 'DIC-TSS'
                        else:
                            nterminal_pblock._annotations.append('upstream start codon spliced out leading to usage of downstream start codon')  # TODO: indicate anchor exon
                            nterminal_pblock.event = 'DIC-splice'
                    else:
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.UP_TIS
                        if caused_by_alt_tss:
                            nterminal_pblock._annotations.append('usage of upstream alternative TIS revealed by alternative TSS')  # TODO: indicate surrounding anchor exons
                            for tblock in nterminal_pblock.transcript_blocks:
                                tblock.flags |= AnnotationFlag.TSS
                        else:
                            nterminal_pblock._annotations.append('usage of upstream alternative TIS revealed by splicing')  # TODO: indicate surrounding anchor exons
                        nterminal_pblock.event = 'upTIS'

                else:
                    # mutually exclusive start codons
                    # TODO: detect if downstream transcript's UTR overlaps upstream transcript's CDS
                    for tblock in nterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.MXIC
                    if alternative_tss:
                        nterminal_pblock._annotations.append('alternative TSS leading to mutually exclusive start codons')
                        for tblock in nterminal_pblock.transcript_blocks:
                            tblock.flags |= AnnotationFlag.TSS
                        nterminal_pblock.event = 'MXIC-TSS'
                    else:
                        nterminal_pblock._annotations.append('5\' UTR splicing leading to mutually exclusive start codons')
                        nterminal_pblock.event = 'MXIC-splice'
        
        # classify C-terminal changes (if any)
        if cterminal_pblock and cterminal_pblock.category is not ProteinAlignCat.MATCH:
            if upstream_cterm_res_aln.category in FRAMESHIFT:
                for tblock in cterminal_pblock.transcript_blocks:
                    tblock.flags |= AnnotationFlag.SIF
                if upstream_cterm_res_aln.anchor.amino_acid is AminoAcid.STOP:
                    cterminal_pblock._annotations.append('splicing-induced frameshift leading to usage of downstream stop codon')  # TODO: indicate location of other stop codon
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.DTC
                    cterminal_pblock.event = 'SIF-DTC'
                else:
                    cterminal_pblock._annotations.append('splicing-induced frameshift leading to usage of upstream stop codon')  # TODO: indicate location of other stop codon
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.UTC
                    cterminal_pblock.event = 'SIF-UTC'
            elif upstream_cterm_res_aln.category is TranscriptAlignCat.DELETION:
                if strand is Strand.PLUS:
                    alt_cterm_exons = upstream_cterm_res_aln.anchor.exons[-1].stop < self.other.orf.exons[-1].start
                elif strand is Strand.MINUS:
                    alt_cterm_exons = upstream_cterm_res_aln.anchor.exons[-1].start > self.other.orf.exons[-1].stop
                if alt_cterm_exons:
                    cterminal_pblock._annotations.append('alternative C-terminal exon')
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.ACTE
                    cterminal_pblock.event = 'ACTE'
                else:
                    cterminal_pblock._annotations.append('upstream stop codon spliced out leading to usage of downstream stop codon')  # TODO: indicate location of other stop codon
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.DTC
                    cterminal_pblock.event = 'DTC-splice'
            elif upstream_cterm_res_aln.category is TranscriptAlignCat.INSERTION:
                exon_extension_introduces_stop = (upstream_cterm_tblock._prev_match_or_frame_tblock and 
                    upstream_cterm_tblock._prev_match_or_frame_tblock[-1].other.codon[0].exon is upstream_cterm_res_aln.other.codon[2].exon)
                if strand is Strand.PLUS:
                    alt_cterm_exons = upstream_cterm_res_aln.other.exons[-1].stop < self.anchor.orf.exons[-1].start
                elif strand is Strand.MINUS:
                    alt_cterm_exons = upstream_cterm_res_aln.other.exons[-1].start > self.anchor.orf.exons[-1].stop
                if exon_extension_introduces_stop:
                    lengthened_exon = upstream_cterm_tblock._prev_match_or_frame_tblock[-1].anchor.codon[0].exon
                    cterminal_pblock._annotations.append(f'upstream stop codon introduced by extension of {lengthened_exon}')
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.EXITC
                    cterminal_pblock.event = 'EXITC'
                elif alt_cterm_exons:
                    cterminal_pblock._annotations.append('alternative C-terminal exon')
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.ACTE
                    cterminal_pblock.event = 'ACTE'
                else:
                    cterminal_pblock._annotations.append('upstream stop codon introduced by splicing')  # TODO: indicate surrounding anchor exons
                    for tblock in cterminal_pblock.transcript_blocks:
                        tblock.flags |= AnnotationFlag.UTC
                    cterminal_pblock.event = 'UTC-splice'
            else:
                cterminal_pblock._annotations.append('complex C-terminal event')
                cterminal_pblock.event = 'complex'
            if self.other.orf.nmd:
                cterminal_pblock._annotations.append('NMD candidate')
                cterminal_pblock.event = 'NMD'

    def get_anchor_residue_alignment(self, anchor_residue: 'Residue', include_index: bool = False):
        if anchor_residue.protein is not self.anchor:
            raise ValueError(f'{anchor_residue} does not belong to {self.anchor}')
        for i, res_aln in enumerate(self):
            if res_aln.anchor is anchor_residue:
                return (i, res_aln) if include_index else res_aln

    def get_other_residue_alignment(self, other_residue: 'Residue', include_index: bool = False):
        if other_residue.protein is not self.other:
            raise ValueError(f'{other_residue} does not belong to {self.other}')
        for i, res_aln in enumerate(self):
            if res_aln.other is other_residue:
                return (i, res_aln) if include_index else res_aln

    def project_feature(self, anchor_feature: 'ProteinFeature'):
        if anchor_feature.protein is not self.anchor:
            raise ValueError(f'{anchor_feature} is not a feature of {self.anchor}')
        feat_aln = FeatureAlignment(self, anchor_feature)
        projected_feature = feat_aln.projected_feature
        return projected_feature, feat_aln
    
    # @property
    # def feature_alignments(self):
    #     return [self.get_feature_alignment(feature) for feature in self.anchor.features]

    # @property
    # def projected_features(self):
    #     return list(filter(None, (feat_aln.projected_feature for feat_aln in self.feature_alignments)))


class FeatureAlignment(AlignmentBlock):
    def __init__(self, parent: 'Alignment', proteinfeature: 'ProteinFeature'):
        first_res = proteinfeature.residues[0]
        last_res = proteinfeature.residues[-1]
        first_res_aln_index, _ = parent.get_anchor_residue_alignment(first_res, include_index=True)
        last_res_aln_index, _ = parent.get_anchor_residue_alignment(last_res, include_index=True)
        super().__init__(parent, first_res_aln_index, last_res_aln_index+1)
        self.proteinfeature = proteinfeature
        self.transcript_blocks = [tblock for tblock in parent.transcript_blocks if tblock.start <= self.end and self.start <= tblock.end]
        self.protein_blocks = [pblock for pblock in parent.protein_blocks if pblock.start <= self.end and self.start <= pblock.end]

    def __repr__(self):
        return f'{self.proteinfeature}>>{self.parent.other}'
    
    @cached_property
    def projected_feature(self):
        if any(res_aln.other for res_aln in self if res_aln.category in MATCH_OR_FRAME):
            residues = self.other_residues
            diff_string = ''.join(str(res_aln.category) for res_aln in self if res_aln.category is not TranscriptAlignCat.DELETION)
            assert len(diff_string) == len(residues)
            proj_feat = ProjectedFeature(
                feature = self.proteinfeature.feature,
                protein = self.parent.other,
                protein_start = residues[0].position,
                protein_stop = residues[-1].position,
                reference = False,
                anchor = self.proteinfeature,
                _differences = run_length_encode(diff_string)
            )
            return proj_feat
        else:
            return None


# class ProjectedFeature:
#     def __init__(self, feature_alignment: 'FeatureAlignment'):
#         self.alignment = feature_alignment

#         self.residues = feature_alignment.other_residues
#         self.altered_residues = [res_aln.other for res_aln in feature_alignment if res_aln.category not in {TranscriptAlignCat.MATCH, TranscriptAlignCat.EDGE_MATCH, TranscriptAlignCat.DELETION}]
#         self.sequence = ''.join(str(res.amino_acid) for res in self.residues)
#         self.protein_start = self.residues[0].position
#         self.protein_stop = self.residues[-1].position

#         self.type = feature_alignment.feature.type
#         self.name = feature_alignment.feature.name
#         self.accession = feature_alignment.feature.accession
#         self.protein_id = feature_alignment.feature.protein_id
#         self.protein = feature_alignment.feature.protein
    
#     def __repr__(self):
#         return f'{self.alignment.feature.protein}>>' + ProteinFeature.__repr__(self)


def rough_alignment(anchor: 'Protein', other: 'Protein', strand: 'Strand') -> List['ResidueAlignment']:
    # helper functions
    def make_gap_res_from_prev_res(prev_res: 'Residue', protein: 'Protein'):
        gap_pos = prev_res.position if prev_res else 0
        upstream_exon = prev_res.codon[0].exon if prev_res else None
        return GapResidue(protein, gap_pos, upstream_exon, None)
    def get_next_res_from_stack(stack: MutableSequence['Residue'], gap_res: 'GapResidue'):
        next_res = stack.popleft()
        if gap_res:
            gap_res.downstream_exon = next_res.codon[-1].exon
        return next_res

    anchor_stack = deque(anchor.residues)
    other_stack = deque(other.residues)
    anchor_res, other_res = None, None
    anchor_gap, other_gap = None, None
    chain = []
    while anchor_stack or other_stack:
        # use residues at top of stack to determine event type
        event_type = TranscriptAlignCat.UNKNOWN
        if anchor_stack and other_stack:
            anchor_current, other_current = anchor_stack[0], other_stack[0]
            if len(anchor_current.exons) == 1 and len(other_current.exons) == 1:
                # for contiguous (nonsplit) codons, can just compare the coords of their middle nucleotides
                coord_diff = anchor_current.codon[1].coordinate - other_current.codon[1].coordinate
                if strand is Strand.MINUS:
                    coord_diff = -coord_diff
                if coord_diff == 0:
                    if anchor_current.amino_acid is other_current.amino_acid:
                        event_type = TranscriptAlignCat.MATCH
                elif coord_diff == -1:
                    event_type = TranscriptAlignCat.FRAME_AHEAD
                elif coord_diff == 1:
                    event_type = TranscriptAlignCat.FRAME_BEHIND
                elif coord_diff < 0:
                    event_type = TranscriptAlignCat.DELETION
                elif coord_diff > 0:
                    event_type = TranscriptAlignCat.INSERTION
            else:
                # for split codons, need to compare the coords of all their nucleotides
                anchor_coords = tuple(nt.coordinate for nt in anchor_current.codon)
                other_coords = tuple(nt.coordinate for nt in other_current.codon)
                overlap = len(set(anchor_coords) & set(other_coords))
                coord_diff = anchor_coords[1] - other_coords[1]
                if strand is Strand.MINUS:
                    coord_diff = -coord_diff
                if overlap < 2:
                    if coord_diff < 0:
                        event_type = TranscriptAlignCat.DELETION
                    elif coord_diff > 0:
                        event_type = TranscriptAlignCat.INSERTION
                elif overlap == 2:
                    if coord_diff < 0:
                        event_type = TranscriptAlignCat.FRAME_AHEAD
                    elif coord_diff > 0:
                        event_type = TranscriptAlignCat.FRAME_BEHIND
                    else:
                        event_type = TranscriptAlignCat.EDGE
                elif overlap == 3 and anchor_current.amino_acid is other_current.amino_acid:
                    event_type = TranscriptAlignCat.MATCH
        elif not other_stack:
            event_type = TranscriptAlignCat.DELETION
        elif not anchor_stack:
            event_type = TranscriptAlignCat.INSERTION

        if event_type is TranscriptAlignCat.DELETION:
            anchor_res = get_next_res_from_stack(anchor_stack, anchor_gap)
            anchor_gap = None
            if not other_gap:
                other_gap = make_gap_res_from_prev_res(other_res, other)
            res_align = ResidueAlignment(anchor_res, other_gap, TranscriptAlignCat.DELETION)
        elif event_type is TranscriptAlignCat.INSERTION:
            other_res = get_next_res_from_stack(other_stack, other_gap)
            other_gap = None
            if not anchor_gap:
                anchor_gap = make_gap_res_from_prev_res(anchor_res, anchor)
            res_align = ResidueAlignment(anchor_gap, other_res, TranscriptAlignCat.INSERTION)
        else:
            anchor_res = get_next_res_from_stack(anchor_stack, anchor_gap)
            other_res = get_next_res_from_stack(other_stack, other_gap)
            anchor_gap, other_gap = None, None
            res_align = ResidueAlignment(anchor_res, other_res, event_type)
        
        assert res_align.anchor.protein is anchor
        assert res_align.other.protein is other
        chain.append(res_align)
    return chain


def refine_alignment(chain: List['ResidueAlignment']) -> None:
    """Performs a second pass on an existing ResidueAlignment chain to identify complex codon alignments."""
    DELETE_INSERT = {TranscriptAlignCat.DELETION, TranscriptAlignCat.INSERTION}
    FRAMESHIFT = {TranscriptAlignCat.FRAME_AHEAD, TranscriptAlignCat.FRAME_BEHIND}
    i = 1
    while i+1 < len(chain):
        curr = chain[i]
        curr_cat = curr.category
        if curr_cat in DELETE_INSERT:
            prev = chain[i-1]
            next = chain[i+1]
            prev_cat = prev.category
            next_cat = next.category
            if len({prev_cat, curr_cat, next_cat}) == 3:
                merge_res = None
                if prev_cat in DELETE_INSERT and next_cat in FRAMESHIFT:
                    merge_res = prev
                elif next_cat in DELETE_INSERT and prev_cat in FRAMESHIFT:
                    merge_res = next
                if merge_res:
                    if curr_cat is TranscriptAlignCat.DELETION:
                        chain[i] = ResidueAlignment(curr.anchor, merge_res.other, TranscriptAlignCat.COMPLEX)
                    elif curr_cat is TranscriptAlignCat.INSERTION:
                        chain[i] = ResidueAlignment(merge_res.anchor, curr.other, TranscriptAlignCat.COMPLEX)
                    if merge_res is prev:
                        chain.pop(i-1)
                        continue  # popping previous item from chain automatically advances i to the next item
                    elif merge_res is next:
                        chain.pop(i+1)
        i += 1


def get_transcript_blocks(aln: Iterable['ResidueAlignment']) -> List['TranscriptAlignmentBlock']:
    tblocks = []
    start = 0
    prev_match_or_frame_tblock = None
    for i, (category, res_alns) in enumerate(groupby(aln, key=attrgetter('category'))):
        tblock_length = len(list(res_alns))
        tblock = TranscriptAlignmentBlock(aln, i, start, start+tblock_length, category)
        if category in MATCH_OR_FRAME:
            for prev_tblock in reversed(tblocks):
                if prev_tblock is prev_match_or_frame_tblock:
                    break
                prev_tblock._next_match_or_frame_tblock = tblock
            prev_match_or_frame_tblock = tblock
        else:
            tblock._prev_match_or_frame_tblock = prev_match_or_frame_tblock
        tblocks.append(tblock)
        start += tblock_length
    return tblocks


def get_protein_blocks(parent: 'Alignment') -> List['ProteinAlignmentBlock']:
    # TODO: account for amino acid sequence
    pblocks = []
    for i, (is_match, tblock_group) in enumerate(groupby(parent.transcript_blocks, key=lambda tblock: tblock.category is TranscriptAlignCat.MATCH)):
        tblock_group = list(tblock_group)
        if is_match:
            pblock_category = ProteinAlignCat.MATCH
        else:
            categories = {tblock.category for tblock in tblock_group if tblock.category is not TranscriptAlignCat.EDGE}
            pblock_category = ProteinAlignCat.SUBSTITUTION
            if len(categories) == 1:
                single_category = list(categories)[0]
                if single_category is TranscriptAlignCat.DELETION:
                    pblock_category = ProteinAlignCat.DELETION
                elif single_category is TranscriptAlignCat.INSERTION:
                    pblock_category = ProteinAlignCat.INSERTION
        pblocks.append(ProteinAlignmentBlock(parent, i, tblock_group, pblock_category))
    return pblocks
    
### helper functions ###

def pairwise_align_protein_sets(setA: Iterable['Protein'], setB: Iterable['Protein']):
    return [Alignment(protA, protB) for protA, protB in product(setA, setB)]


def export_annotated_pblocks_to_tsv(output_path, pblocks: Iterable['ProteinAlignmentBlock']):
    with open(output_path, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=PBLOCK_FIELDS, delimiter='\t', quotechar='"')
        writer.writeheader()
        for pblock in pblocks:
            if pblock.annotation:
                writer.writerow(pblock.to_dict())
