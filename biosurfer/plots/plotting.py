# functions to create different visualizations of isoforms/clones/domains/muts
from copy import copy
from dataclasses import dataclass
from itertools import chain, groupby, islice
from operator import attrgetter, sub
from typing import (TYPE_CHECKING, Any, Callable, Collection, Dict, Iterable,
                    List, Literal, Optional, Set, Tuple, Union)
from warnings import filterwarnings, warn
from Bio import Align

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from biosurfer.core.alignments import (FRAMESHIFT, ProjectedFeature,
                                       Alignment)
from biosurfer.core.constants import (AminoAcid, FeatureType, SequenceAlignmentCategory,
                                      Strand, CodonAlignmentCategory)
from biosurfer.core.helpers import ExceptionLogger, Interval, IntervalTree, get_interval_overlap_graph
from biosurfer.core.models.biomolecules import (GencodeTranscript,
                                                PacBioTranscript, Transcript)
from brokenaxes import BrokenAxes
from graph_tool import Graph
from graph_tool.topology import sequential_vertex_coloring
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

if TYPE_CHECKING:
    from biosurfer.core.alignments import ProteinAlignmentBlock
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

StartStop = Tuple[int, int]

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

# colors for different transcript types
TRANSCRIPT_COLORS = {
    GencodeTranscript: ('#343553', '#5D5E7C'),
    PacBioTranscript: ('#61374D', '#91677D')
}

# alpha values for different absolute reading frames
ABS_FRAME_ALPHA = {0: 1.0, 1: 0.45, 2: 0.15}

# hatching styles for different relative frameshifts
REL_FRAME_STYLE = {
    CodonAlignmentCategory.FRAME_AHEAD: '////',
    CodonAlignmentCategory.FRAME_BEHIND: 'xxxx'
}

PBLOCK_COLORS = {
    SequenceAlignmentCategory.DELETION: '#FF0082',
    SequenceAlignmentCategory.INSERTION: '#05E0FF',
    SequenceAlignmentCategory.SUBSTITUTION: '#FFD700'
}

FEATURE_COLORS = {
    'MobiDB': '#AAAAAA'
}

@dataclass
class IsoformPlotOptions:
    """Bundles various options for adjusting plots made by IsoformPlot."""
    intron_spacing: int = 30  # number of bases to show in each intron
    track_spacing: float = 1.5  # ratio of space between tracks to max track width
    subtle_splicing_threshold: int = 20  # maximum difference (in bases) between exon boundaries to display subtle splicing

    @property
    def max_track_width(self) -> float:
        return 1/(self.track_spacing + 1)
    
    @max_track_width.setter
    def max_track_width(self, width: float):
        self.track_spacing = (1 - width)/width


TableColumn = Callable[[Transcript], str]
class IsoformPlot:
    """Encapsulates methods for drawing one or more isoforms aligned to the same genomic x-axis."""
    def __init__(self, transcripts: Iterable['Transcript'], columns: Dict[str, TableColumn] = None, **kwargs):
        self.transcripts: List['Transcript'] = list(transcripts)  # list of orf objects to be drawn
        gene = {tx.gene for tx in filter(None, self.transcripts)}
        if len(gene) > 1:
            raise ValueError(f'Found isoforms from multiple genes: {", ".join(g.name for g in gene)}')
        strand = {tx.strand for tx in filter(None, self.transcripts)}
        if len(strand) > 1:
            raise ValueError("Can't plot isoforms from different strands")
        self.strand: Strand = list(strand)[0]

        self.fig: Optional['Figure'] = None
        self._bax: Optional['BrokenAxes'] = None
        self._columns: Dict[str, TableColumn] = columns if columns else {'': lambda x: ''}
        self.opts = IsoformPlotOptions(**kwargs)
        self.reset_xlims()

        # keep track of artists for legend
        self._handles = dict()


    # Internally, IsoformPlot stores _subaxes, which maps each genomic region to the subaxes that plots the region's features.
    # The xlims property provides a simple interface to allow users to control which genomic regions are plotted.
    @property
    def xlims(self) -> Tuple[StartStop]:
        """Coordinates of the genomic regions to be plotted, as a tuple of (start, end) tuples."""
        return self._xlims

    @xlims.setter
    def xlims(self, xlims: Iterable[StartStop]):
        xregions = IntervalTree.from_tuples((min(start, stop), max(start, stop)+1) for start, stop in xlims)  # condense xlims into single IntervalTree object
        xregions.merge_equals()
        xregions.merge_overlaps()
        xregions.merge_neighbors()
        xregions = sorted(xregions.all_intervals)
        if self.strand is Strand.MINUS:
            xregions.reverse()
        self._subaxes = IntervalTree.from_tuples((start, end, i) for i, (start, end, _) in enumerate(xregions))
        self._xlims = tuple((start, end-1) if self.strand is Strand.PLUS else (end-1, start) for start, end , _ in xregions)

    def reset_xlims(self):
        """Set xlims automatically based on exons in isoforms."""
        space = self.opts.intron_spacing
        self.xlims = tuple((exon.start - space, exon.stop + space) for tx in filter(None, self.transcripts) for exon in tx.exons)

    # This method speeds up plotting by allowing IsoformPlot to add artists only to the subaxes where they are needed.
    def _get_subaxes(self, xcoords: Union[int, StartStop]) -> Tuple['Axes']:
        """For a specific coordinate or range of coordinates, retrieve corresponding subaxes."""
        if isinstance(xcoords, tuple):
            if xcoords[0] > xcoords[1]:
                xcoords = (xcoords[1], xcoords[0])
            xcoords = slice(*xcoords)
        subax_ids = [interval[-1] for interval in self._subaxes[xcoords]]
        if not subax_ids:
            raise ValueError(f"{xcoords} is not within plot's xlims")
        return tuple(self._bax.axs[id] for id in subax_ids)

    def draw_point(self, track: int, pos: int,
                    y_offset: float = 0.0,
                    height: Optional[float] = None,
                    type='line', marker='.', linewidth=1, **kwargs):
        """Draw a feature at a specific point. Appearance types are 'line' and 'lollipop'."""
        # TODO: make type an enum?
        if type == 'line':
            if height is None:
                height = self.opts.max_track_width
            center = track + y_offset
            artist = mlines.Line2D(
                xdata = (pos, pos),
                ydata = (center - height/2, center + height/2),
                linewidth = linewidth,
                **kwargs
            )
        elif type == 'lollipop':
            if height is None:
                height = 0.3*self.opts.max_track_width
            artist = mlines.Line2D(
                xdata = (pos, pos),
                ydata = (-0.25 - height, -0.25),
                linewidth = linewidth,
                marker = marker,
                markevery = 2,
                **kwargs
            )
        else:
            raise ValueError(f'Point type "{type}" is not defined')
        
        try:
            subaxes = self._get_subaxes(pos)[0]
        except ValueError as e:
            warn(str(e))
        else:
            subaxes.add_artist(artist)
        return artist

    def draw_region(self, track: int, start: int, stop: int,
                    y_offset: Optional[float] = None,
                    height: Optional[float] = None,
                    type='rect', **kwargs):
        """Draw a feature that spans a region. Appearance types are rectangle and line."""
        if start == stop:
            return
        # TODO: make type an enum?
        if type == 'rect':
            if height is None:
                height = self.opts.max_track_width
            if y_offset is None:
                y_offset = -0.5*height
            artist = mpatches.Rectangle(
                xy = (start, track + y_offset),
                width = stop - start,
                height = height,
                **kwargs
            )
        elif type == 'line':
            if y_offset is None:
                y_offset = 0
            artist = mlines.Line2D(
                xdata = (start, stop),
                ydata = (track + y_offset, track + y_offset),
                **kwargs
            )
        else:
            raise ValueError(f'Region type "{type}" is not defined')

        subaxes = self._get_subaxes((start, stop))
        for ax in subaxes:
            ax.add_artist(copy(artist))
        return artist
    
    def draw_background_rect(self, start: int, stop: int,
                            track_first: int = None, track_last: int = None,
                            padding: float = None, **kwargs):
        """Draw a rectangle in the background of the plot."""
        if start == stop:
            return
        if track_first is None:
            track_first = 0
        if track_last is None:
            track_last = len(self.transcripts) - 1
        if padding is None:
            padding = self.opts.max_track_width
        top = track_first - padding
        bottom = track_last + padding
        artist = mpatches.Rectangle(
            xy = (start, top),
            width = stop - start,
            height = bottom - top,
            zorder = 0.5,
            **kwargs
        )

        subaxes = self._get_subaxes((start, stop))
        for ax in subaxes:
            ax.add_artist(copy(artist))
        return artist

    def draw_text(self, x: int, y: float, text: str, **kwargs):
        """Draw text at a specific location. x-coordinate is genomic, y-coordinate is w/ respect to tracks (0-indexed).
        Ex: x=20000, y=2 will center text on track 2 at position 20,000."""
        # TODO: make this use Axes.annotate instead
        # we can't know how much horizontal space text will take up ahead of time
        # so text is plotted using BrokenAxes' big_ax, since it spans the entire x-axis
        big_ax = self._bax.big_ax
        try:
            subaxes = self._get_subaxes(x)[0]  # grab coord transform from correct subaxes
        except ValueError as e:
            warn(str(e))
        else:
            big_ax.text(x, y, text, transform=subaxes.transData, **kwargs)

    def draw_legend(self, only_labels: Optional[Iterable[str]] = None, except_labels: Optional[Iterable[str]] = None, **kwargs):
        if only_labels and except_labels:
            raise ValueError('Cannot set both "only_labels" and "except_labels"')
        elif only_labels:
            labels = [label for label in self._handles if label in only_labels]
        elif except_labels:
            labels = [label for label in self._handles if label not in except_labels]
        else:
            labels = list(self._handles.keys())
        handles = [self._handles[label] for label in labels]
        self.fig.legend(
            handles = handles,
            labels = labels,
            # ncol = 1,
            # loc = 'center left',
            # mode = 'expand',
            # bbox_to_anchor = (1.05, 0.5),
            **kwargs
        )

    def draw_isoform(self, tx: 'Transcript', track: int):
        """Plot a single isoform in the given track."""
        start, stop = tx.start, tx.stop
        align_start, align_stop = 'right', 'left'
        if self.strand is Strand.MINUS:
            align_start, align_stop = align_stop, align_start
        
        # plot intron line
        self.draw_region(
            track,
            start = start,
            stop = stop,
            type = 'line',
            linewidth = 1.5,
            color = 'gray',
            zorder = 1.5
        )

        # plot exons
        utr_kwargs = {
            'type': 'rect',
            'edgecolor': 'k',
            'facecolor': TRANSCRIPT_COLORS[type(tx)][1],
            'height': 0.5*self.opts.max_track_width,
            'zorder': 1.5
        }
        cds_kwargs = {
            'type': 'rect',
            'edgecolor': 'k',
            'facecolor': TRANSCRIPT_COLORS[type(tx)][0],
            'zorder': 1.5
        }
        if tx.orfs:
            orf = tx.primary_orf
            if orf.utr5:
                for exon in orf.utr5.exons:
                    if self.strand is Strand.PLUS:
                        start = exon.start
                        stop = min(exon.stop, orf.start)
                    elif self.strand is Strand.MINUS:
                        start = max(exon.start, orf.stop)
                        stop = exon.stop
                    self.draw_region(track, start=start, stop=stop, **utr_kwargs)
            for exon in orf.exons:
                start = max(exon.start, orf.start)
                stop = min(exon.stop, orf.stop)
                self.draw_region(track, start=start, stop=stop, **cds_kwargs)
            if orf.utr3:
                for exon in orf.utr3.exons:
                    if self.strand is Strand.PLUS:
                        start = max(exon.start, orf.stop)
                        stop = exon.stop
                    elif self.strand is Strand.MINUS:
                        start = exon.start
                        stop = min(exon.stop, orf.start)
                    self.draw_region(track, start=start, stop=stop, **utr_kwargs)
        else:
            for exon in tx.exons:
                self.draw_region(track, start=exon.start, stop=exon.stop, **utr_kwargs)
        
        for exon in tx.exons:
            # label every 5th exon in anchor isoform for easier navigation
            if track == 0 and exon.position % 5 == 0:
                self.draw_text((exon.start + exon.stop)//2, track - self.opts.max_track_width, f'E{exon.position}', ha='center', va='baseline')            
                        
            # add subtle splice (delta) amounts, if option turned on
            # first, make sure the exon contains a (coding) cds object
            # if exon.cds:
            #     # TODO: pull subtle splice detection code out into this method?
            #     delta_start, delta_end = retrieve_subtle_splice_amounts(exon.cds)
            #     if delta_start:
            #         self.draw_text(exon.start, track-0.1, delta_start, va='bottom', ha=align_start, size='x-small')
            #     if delta_end:
            #         self.draw_text(exon.stop, track-0.1, delta_end, va='bottom', ha=align_stop, size='x-small')
        
        for orf in tx.orfs:
            first_res = orf.protein.residues[0]
            last_res = orf.protein.residues[-1]
            if first_res.amino_acid is AminoAcid.MET:
                start_codon = first_res.codon[0].coordinate
                self.draw_point(track, start_codon, type='line', color='lime')
            if last_res.amino_acid is AminoAcid.STOP:
                stop_codon = last_res.codon[2].coordinate
                self.draw_point(track, stop_codon, type='line', color='red')

        if hasattr(tx, 'start_nf') and tx.start_nf:
            self.draw_text(tx.start if self.strand is Strand.PLUS else tx.stop, track, '! ', ha='right', va='center', weight='bold', color='r')
        if hasattr(tx, 'end_nf') and tx.end_nf:
            self.draw_text(tx.stop if self.strand is Strand.PLUS else tx.start, track, ' !', ha='left', va='center', weight='bold', color='r')

    def draw_all_isoforms(self, subplot_spec = None):
        """Plot all isoforms."""
        R = len(self.transcripts)
        C = len(self._columns)
        self.fig = plt.figure()
        self._bax = BrokenAxes(fig=self.fig, xlims=self.xlims, ylims=((R-0.5, -0.5),), wspace=0, d=0.008, subplot_spec=subplot_spec)
        self._handles['intron'] = mlines.Line2D([], [], linewidth=1.5, color='gray')
        self._handles['exon (GENCODE)'] = mpatches.Patch(facecolor=TRANSCRIPT_COLORS[GencodeTranscript][0], edgecolor='k')
        self._handles['exon (PacBio)'] = mpatches.Patch(facecolor=TRANSCRIPT_COLORS[PacBioTranscript][0], edgecolor='k')
        self._handles['start codon'] = mlines.Line2D([], [], linestyle='None', color='lime', marker='|', markersize=10, markeredgewidth=1)
        self._handles['stop codon'] = mlines.Line2D([], [], linestyle='None', color='red', marker='|', markersize=10, markeredgewidth=1)


        # process orfs to get ready for plotting
        # find_and_set_subtle_splicing_status(self.transcripts, self.opts.subtle_splicing_threshold)
        
        for i, tx in enumerate(self.transcripts):
            with ExceptionLogger(f'Error plotting {tx}'):
                if tx:
                    self.draw_isoform(tx, i)
        
        # plot genomic region label
        # gene = self.transcripts[0].gene
        # start, end = self.xlims[0][0], self.xlims[-1][1]
        # self._bax.set_title(f'{gene.chromosome}({self.strand}):{start}-{end}')
        
        # hide y axis spine
        left_subaxes = self._bax.axs[0]
        left_subaxes.spines['left'].set_visible(False)
        left_subaxes.set_yticks([])
        
        # plot table
        # https://stackoverflow.com/a/57169705
        table = self._bax.big_ax.table(
            rowLabels = [getattr(tx, 'name', '') for tx in self.transcripts],
            colLabels = list(self._columns.keys()),
            cellText = [[f(tx) if tx else '' for f in self._columns.values()] for tx in self.transcripts],
            cellLoc = 'center',
            edges = 'open',
            bbox = (-0.1*C, 0.0, 0.1*C, (R+1)/R)
        )
        # table.auto_set_font_size(False)
        # table.set_fontsize(10)
        
        # rotate x axis tick labels for better readability
        for subaxes in self._bax.axs:
            subaxes.xaxis.set_major_formatter('{x:.0f}')
            for label in subaxes.get_xticklabels():
                label.set_va('top')
                label.set_rotation(90)
                label.set_size(8)
    
    def draw_frameshifts(self, anchor: Optional['Transcript'] = None, hatch_color='white'):
        """Plot relative frameshifts on all isoforms. Uses first isoform as the anchor by default."""
        self._handles['frame +1'] = mpatches.Patch(facecolor='k', edgecolor='w', hatch=REL_FRAME_STYLE[CodonAlignmentCategory.FRAME_AHEAD])
        self._handles['frame -1'] = mpatches.Patch(facecolor='k', edgecolor='w', hatch=REL_FRAME_STYLE[CodonAlignmentCategory.FRAME_BEHIND])
        
        if anchor is None:
            anchor = next(filter(None, self.transcripts))
        if not anchor or not anchor.protein:
            warn(
                'Cannot draw frameshifts without an anchor ORF'
            )
            return
        for i, other in enumerate(self.transcripts):
            if not other or not other.protein or other is anchor:
                continue
            aln = Alignment(anchor.protein, other.protein)
            for (category, exons), block in groupby(aln, key=attrgetter('category', 'other.exons')):
                if category in FRAMESHIFT:
                    if len(exons) > 1:
                        continue
                    block = list(block)
                    start = block[0].other.codon[0].coordinate
                    stop = block[-1].other.codon[2].coordinate
                    self.draw_region(
                        track = i,
                        start = start,
                        stop = stop,
                        facecolor = 'none',
                        edgecolor = hatch_color,
                        linewidth = 0.0,
                        zorder = 1.9,
                        hatch = REL_FRAME_STYLE[category]
                    )
    
    def draw_protein_block(self, pblock: 'ProteinAlignmentBlock', alpha: float = 0.5):
        if 'deletion' not in self._handles:
            self._handles['deletion'] = mpatches.Patch(facecolor=PBLOCK_COLORS[SequenceAlignmentCategory.DELETION])
        if 'insertion' not in self._handles:
            self._handles['insertion'] = mpatches.Patch(facecolor=PBLOCK_COLORS[SequenceAlignmentCategory.INSERTION])
        if 'substitution' not in self._handles:
            self._handles['substitution'] = mpatches.Patch(facecolor=PBLOCK_COLORS[SequenceAlignmentCategory.SUBSTITUTION])

        if pblock.category is SequenceAlignmentCategory.DELETION:
            other_start = pblock.anchor_residues[0].codon[1].coordinate
            other_stop = pblock.anchor_residues[-1].codon[1].coordinate
        else:
            other_start = pblock.other_residues[0].codon[1].coordinate
            other_stop = pblock.other_residues[-1].codon[1].coordinate
        if pblock.category is SequenceAlignmentCategory.INSERTION:
            anchor_start = pblock.other_residues[0].codon[1].coordinate
            anchor_stop = pblock.other_residues[-1].codon[1].coordinate
        else:
            anchor_start = pblock.anchor_residues[0].codon[1].coordinate
            anchor_stop = pblock.anchor_residues[-1].codon[1].coordinate
        self.draw_region(
            self.transcripts.index(pblock.parent.anchor.transcript),
            start = anchor_start,
            stop = anchor_stop,
            y_offset = -0.9*self.opts.max_track_width,
            height = 0.4*self.opts.max_track_width,
            edgecolor = 'none',
            facecolor = PBLOCK_COLORS[pblock.category],
            alpha = alpha
        )
        self.draw_region(
            self.transcripts.index(pblock.parent.other.transcript),
            start = other_start,
            stop = other_stop,
            y_offset = 0.5*self.opts.max_track_width,
            height = 0.4*self.opts.max_track_width,
            edgecolor = 'none',
            facecolor = PBLOCK_COLORS[pblock.category],
            alpha = alpha
        )
    
    def draw_features(self):
        h = self.opts.max_track_width
        feature_names = sorted({feature.name for tx in filter(None, self.transcripts) if tx.protein for feature in tx.protein.features if feature.type is not FeatureType.IDR})
        cmap = sns.color_palette('pastel', len(feature_names))
        colors = dict(zip(feature_names, cmap))
        colors.update(FEATURE_COLORS)
        self._handles.update({name: mpatches.Patch(facecolor=color) for name, color in colors.items()})
        for track, tx in enumerate(self.transcripts):
            if not tx or not tx.protein:
                continue
            features = tx.protein.features
            if not features:
                continue
            subtracks, n_subtracks = generate_subtracks(
                ((feature.protein_start, feature.protein_stop) for feature in features),
                (feature.name for feature in features)
            )
            for feature in features:
                subtrack = subtracks[feature.name]
                color = colors[feature.name]
                if feature.reference:
                    subfeatures = groupby(feature.residues, key=lambda res: (False, res.primary_exon))
                    # n_subtracks_temp = n_subtracks
                else:
                    subfeatures = groupby(feature.residues, key=lambda res: (res in feature.altered_residues, res.primary_exon))
                    # n_subtracks_temp = 2*n_subtracks
                for (altered, _), subfeature in subfeatures:
                    subfeature = list(subfeature)
                    start = subfeature[0].codon[1].coordinate
                    stop = subfeature[-1].codon[1].coordinate
                    self.draw_region(
                        track,
                        start = start,
                        stop = stop,
                        y_offset = (-0.5 + subtrack/n_subtracks)*h,
                        height = h/n_subtracks,
                        edgecolor = 'none',
                        facecolor = color,
                        alpha = 0.5 if altered else 1.0,
                        zorder = 1.8,
                        label = feature.name
                    )
                # draw box behind entire feature
                self.draw_region(
                    track,
                    start = feature.residues[0].codon[1].coordinate,
                    stop = feature.residues[-1].codon[1].coordinate,
                    y_offset = (-0.5 + subtrack/n_subtracks)*h,
                    height = h/n_subtracks,
                    edgecolor = 'none',
                    facecolor = color,
                    alpha = 0.5,
                    zorder = 1.4
                )


def generate_subtracks(intervals: Iterable[Tuple[int, int]], labels: Iterable):
    # inspired by https://stackoverflow.com/a/19088519
    # build graph of labels where labels are adjacent if their intervals overlap
    g, vertex_labels = get_interval_overlap_graph(intervals, labels)
    # find vertex coloring of graph
    # all labels w/ same color can be put into same subtrack
    coloring = sequential_vertex_coloring(g)
    label_to_subtrack = dict(zip(vertex_labels, coloring))
    subtracks = max(label_to_subtrack.values(), default=0) + 1
    return label_to_subtrack, subtracks
