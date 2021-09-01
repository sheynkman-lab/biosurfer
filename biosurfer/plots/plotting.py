# functions to create different visualizations of isoforms/clones/domains/muts

from copy import copy
from dataclasses import dataclass
from itertools import groupby
from operator import attrgetter, sub
from typing import (TYPE_CHECKING, Collection, Dict, Iterable, List, Literal,
                    Optional, Set, Tuple, Union)

import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from biosurfer.core.alignments import TranscriptBasedAlignment
from biosurfer.core.constants import TranscriptLevelAlignmentCategory
from biosurfer.core.helpers import Interval, IntervalTree
from biosurfer.core.models import (ORF, Gene, Junction, Protein, Strand,
                                   Transcript)
from brokenaxes import BrokenAxes

if TYPE_CHECKING:
    from matplotlib.axes import Axes
    from matplotlib.figure import Figure

StartStop = Tuple[int, int]

# alpha values for different absolute reading frames
ABS_FRAME_ALPHA = {0: 1.0, 1: 0.45, 2: 0.15}

# hatching styles for different relative frameshifts
REL_FRAME_STYLE = {
    TranscriptLevelAlignmentCategory.FRAME_AHEAD: '////',
    TranscriptLevelAlignmentCategory.FRAME_BEHIND: 'xxxx'
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


class IsoformPlot:
    """Encapsulates methods for drawing one or more isoforms aligned to the same genomic x-axis."""
    def __init__(self, transcripts: Iterable['Transcript'], **kwargs):
        self.transcripts: List['Transcript'] = list(transcripts)  # list of orf objects to be drawn
        gene = {tx.gene for tx in self.transcripts}
        if len(gene) > 1:
            raise ValueError(f'Found isoforms from multiple genes: {", ".join(g.name for g in gene)}')
        strand = {tx.strand for tx in self.transcripts}
        if len(strand) > 1:
            raise ValueError("Can't plot isoforms from different strands")
        self.strand: Strand = list(strand)[0]

        self._bax: Optional['BrokenAxes'] = None
        self.opts = IsoformPlotOptions(**kwargs)
        self.reset_xlims()
    
    @property
    def fig(self) -> Optional['Figure']:
        return self._bax.fig if self._bax else None

    # Internally, IsoformPlot stores _subaxes, which maps each genomic region to the subaxes that plots the region's features.
    # The xlims property provides a simple interface to allow users to control which genomic regions are plotted.
    @property
    def xlims(self) -> Tuple[StartStop]:
        """Coordinates of the genomic regions to be plotted, as a tuple of (start, end) tuples."""
        return self._xlims

    @xlims.setter
    def xlims(self, xlims: Iterable[StartStop]):
        xregions = IntervalTree.from_tuples((start, stop+1) for start, stop in xlims)  # condense xlims into single IntervalTree object
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
        self.xlims = tuple((exon.start - space, exon.stop + space) for tx in self.transcripts for exon in tx.exons)

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
        
        subaxes = self._get_subaxes(pos)[0]
        subaxes.add_artist(artist)

    def draw_region(self, track: int, start: int, stop: int,
                    y_offset: Optional[float] = None,
                    height: Optional[float] = None,
                    type='rect', **kwargs):
        """Draw a feature that spans a region. Appearance types are rectangle and line."""
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
    
    def draw_background_rect(self, start: int, stop: int,
                            track_first: int = None, track_last: int = None,
                            padding: float = None, **kwargs):
        """Draw a rectangle in the background of the plot."""
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

    # TODO: implement draw_track_label
    def draw_track_label(self):
        pass

    def draw_text(self, x: int, y: float, text: str, **kwargs):
        """Draw text at a specific location. x-coordinate is genomic, y-coordinate is w/ respect to tracks (0-indexed).
        Ex: x=20000, y=2 will center text on track 2 at position 20,000."""
        # TODO: make this use Axes.annotate instead
        # we can't know how much horizontal space text will take up ahead of time
        # so text is plotted using BrokenAxes' big_ax, since it spans the entire x-axis
        big_ax = self._bax.big_ax
        subaxes = self._get_subaxes(x)[0]  # grab coord transform from correct subaxes
        big_ax.text(x, y, text, transform=subaxes.transData, **kwargs)

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
            'facecolor': 'lightsteelblue',
            'height': 0.5*self.opts.max_track_width,
            'zorder': 1.5
        }
        cds_kwargs = {
            'type': 'rect',
            'edgecolor': 'k',
            'facecolor': 'steelblue',
            'zorder': 1.5
        }
        orf = tx.orfs[0]
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
            start_codon = orf.protein.residues[0].codon[0].coordinate
            stop_codon = orf.protein.residues[-1].codon[2].coordinate
            self.draw_point(track, start_codon, type='line', color='lime')
            self.draw_point(track, stop_codon, type='line', color='red')

        if hasattr(tx, 'start_nf') and tx.start_nf:
            self.draw_text(tx.start if self.strand is Strand.PLUS else tx.stop, track, '! ', ha='right', va='center', weight='bold', color='r')
        if hasattr(tx, 'end_nf') and tx.end_nf:
            self.draw_text(tx.stop if self.strand is Strand.PLUS else tx.start, track, ' !', ha='left', va='center', weight='bold', color='r')

    def draw_all_isoforms(self):
        """Plot all isoforms."""
        
        self._bax = BrokenAxes(xlims=self.xlims, ylims=((len(self.transcripts), -2*self.opts.max_track_width),), wspace=0, d=0.008)

        # process orfs to get ready for plotting
        # find_and_set_subtle_splicing_status(self.transcripts, self.opts.subtle_splicing_threshold)
        
        for i, tx in enumerate(self.transcripts):
            self.draw_isoform(tx, i)
        
        # set title
        gene = self.transcripts[0].gene
        start, end = self.xlims[0][0], self.xlims[-1][1]
        self._bax.set_title(f'{gene.chromosome}({self.strand}):{start}-{end}')
        
        # hide y axis spine and replace tick labels with ORF ids
        left_subaxes = self._bax.axs[0]
        left_subaxes.spines['left'].set_visible(False)
        left_subaxes.set_yticks(list(range(len(self.transcripts))))
        left_subaxes.set_yticklabels([tx.name for tx in self.transcripts])
        
        # rotate x axis tick labels for better readability
        for subaxes in self._bax.axs:
            subaxes.xaxis.set_major_formatter('{x:.0f}')
            for label in subaxes.get_xticklabels():
                label.set_va('top')
                label.set_rotation(90)
                label.set_size(8)
    
    def draw_frameshifts(self, hatch_color = 'white'):
        """Plot relative frameshifts on all isoforms, using the first isoform as the anchor."""
        FRAMESHIFT = {TranscriptLevelAlignmentCategory.FRAME_AHEAD, TranscriptLevelAlignmentCategory.FRAME_BEHIND}
        anchor_tx = self.transcripts[0]
        anchor = anchor_tx.orfs[0].protein
        for i, other_tx in enumerate(self.transcripts[1:], start=1):
            if not other_tx.orfs:
                continue
            other = other_tx.orfs[0].protein
            aln = TranscriptBasedAlignment(anchor, other)
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
                        zorder = 1.5,
                        hatch = REL_FRAME_STYLE[category]
                    )
