# functions to create different visualizations of isoforms/clones/domains/muts
# deleted all mut-related code (see orig. iso-image script to retrieve)

import os
import math
import operator
from itertools import groupby
from .helpers import IntegerInterval as Interval
from portion import IntervalDict

# writing isoimage to write-out to matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
from brokenaxes import BrokenAxes
from copy import copy
from dataclasses import dataclass

from typing import TYPE_CHECKING, Optional, Union, Literal, List, Tuple, Dict, Set, Iterable, Collection
if TYPE_CHECKING:
    from .isoclass import Gene, ORF
    from matplotlib.axes import Axes
Pair = Tuple[int, int]

# alpha values for different absolute reading frames
ABS_FRAME_ALPHA = {0: 1.0, 1: 0.45, 2: 0.15}


class IsoformPlot:
    """Encapsulates methods for drawing one or more isoforms aligned to the same genomic x-axis."""
    def __init__(self, orfs_to_plot: Iterable['ORF'], **kwargs):
        self.orfs: List['ORF'] = list(orfs_to_plot)  # list of orf objects to be drawn
        self.opts = IsoformPlotOptions(**kwargs)

        self.reset_xlims()
    
    # Internally, IsoformPlot stores _subaxes, which maps each genomic region to the subaxes that plots the region's features.
    # The xlims property provides a simple interface to allow users to control which genomic regions are plotted.
    @property
    def xlims(self) -> Tuple[Pair]:
        """Coordinates of the genomic regions to be plotted, as a tuple of (start, end) tuples."""
        return Interval(self._subaxes.domain()).to_tuples()

    @xlims.setter
    def xlims(self, xlims: Collection[Pair]):
        xregions = Interval.from_tuples(*xlims)  # condense xlims into single Interval object
        self._subaxes: 'IntervalDict[Interval, int]' = IntervalDict({region: i for i, region in enumerate(xregions)})

    def reset_xlims(self):
        """Set xlims automatically based on exons in isoforms."""
        space = self.opts.intron_spacing
        self.xlims = tuple((exon.start - space, exon.end + space) for orf in self.orfs for exon in orf.exons)

    def _get_subaxes(self, xcoords: Union[int, Pair]) -> Tuple['Axes']:
        """For a specific coordinate or range of coordinates, retrieve corresponding subaxes."""
        if isinstance(xcoords, int):
            subax_ids = [self._subaxes[xcoords]]
        elif isinstance(xcoords, tuple):
            subax_ids = list(self._subaxes[Interval.from_tuples(xcoords)].values())
        else:
            raise TypeError('xcoords must be an int or a tuple of 2 ints!')
        return tuple(self._bax.axs[i] for i in subax_ids)

    def draw_region(self, track: int, start: int, end: int,
                    y_offset: Optional[float] = None,
                    height: Optional[float] = None,
                    type='rect', **kwargs):
        """Draw a feature that spans a region. Default type is rectangle."""
        if type == 'rect':
            if y_offset is None:
                y_offset = -0.25
            if height is None:
                height = 0.5
            artist = mpatches.Rectangle(
                xy = (start, track + y_offset),
                width = end - start,
                height = height,
                **kwargs
            )
        elif type == 'line':
            if y_offset is None:
                y_offset = 0
            artist = mlines.Line2D(
                xdata = (start, end),
                ydata = (track + y_offset, track + y_offset),
                **kwargs
            )
        else:
            raise ValueError(f'region type {type} is not defined')

        subaxes = self._get_subaxes((start, end))
        for ax in subaxes:
            ax.add_artist(copy(artist))

    # TODO: implement draw_track_label
    def draw_track_label(self):
        pass

    def draw_orf(self, orf: 'ORF', track: int):
        """Plot a single orf in the given track."""

        # plot orf name
        label = retrieve_orf_name(orf)
        # plot intron line
        self.draw_region(
            track,
            start = orf.first.coord,
            end = orf.last.coord,
            type = 'line',
            linewidth = 1.5,
            color = 'k',
            linestyle = '--',
            dashes = (1, 1.2),
            zorder = 1.5
        )

        # render thin exon blocks
        for exon in orf.exons:
            col = get_orf_color(orf)
            alpha_val = ABS_FRAME_ALPHA[exon.abs_frm] if self.opts.show_abs_frame else 1.0
            self.draw_region(
                track,
                start = exon.start,
                end = exon.end,
                type = 'rect',
                edgecolor = 'k',
                facecolor = col,
                zorder = 1.5
            )
            
            # # add subtle splice (delta) amounts, if option turned on
            # # first, make sure the exon contains a (coding) cds object
            # if exon.cds:
            #     delta_start, delta_end = retrieve_subtle_splice_amounts(exon.cds)
            #     # if delta_start:
            #         # self.ax.text(x, y+0.2, delta_start, va='bottom', ha='left', size='x-small')
            #     # if delta_end:
            #         # self.ax.text(x+blength, y+0.2, delta_end, va='bottom', ha='right', size='x-small')

    def draw(self):
        """Plot all orfs."""
        
        self._bax = BrokenAxes(xlims=self.xlims, ylims=((len(self.orfs), -0.5),), wspace=0)

        # process orfs to get ready for plotting
        gen_obj = grab_gen_objs_from_orfs(self.orfs) # all gen_objs into a set
        verify_only_one_gen_obj_represented(gen_obj)
        verify_all_orfs_of_same_strand(self.orfs)
        # compress_introns_and_set_relative_orf_exon_and_cds_coords(gen_obj, self.orfs, self.opts.intron_spacing)
        find_and_set_subtle_splicing_status(self.orfs, self.opts.subtle_threshold)
        
        repr_orf = get_repr_orf(self.orfs)

        for i, orf in enumerate(self.orfs):
            self.draw_orf(orf, i)
        
        # hide y axis spine and replace tick labels with ORF ids
        left_subaxes = self._bax.axs[0]
        left_subaxes.spines['left'].set_visible(False)
        left_subaxes.set_yticks(list(range(len(self.orfs))))
        left_subaxes.set_yticklabels([orf.name for orf in self.orfs])
        
        # rotate x axis tick labels for better readability
        for subaxes in self._bax.axs:
            labels = [int(label) for label in subaxes.get_xticks()]
            subaxes.tick_params(axis='x', labelsize=8)
            subaxes.set_xticklabels(labels, rotation=90, ha='center', va='top')


PlotMode = Literal['all', 'cds']
@dataclass
class IsoformPlotOptions:
    """Bundles various options for adjusting plots made by IsoformPlot."""
    compress: float = 300.0
    height: float = 0.1
    intron_spacing: int = 30  # number of bases to show in each intron
    label_x: float = -0.1
    mode: PlotMode = 'all'
    show_abs_frame: bool = False
    spacing: float = 0.5
    subtle_threshold: int = 20


### older methods ###


def get_union(intervals: Collection[Interval]) -> Tuple[Interval]:
    """Takes collection of (inclusive) intervals and returns their (possibly disjoint) union."""

    def union(a: Interval, b: Interval) -> Interval:
        # only works if a and b are adjacent
        return min(a[0], b[0]), max(a[1], b[1])
    
    def adjacent(a: Interval, b: Interval) -> bool:
        # (1, 3) and (2, 4) are adjacent
        # (1, 3) and (3, 5) are adjacent
        # (1, 3) and (4, 6) are adjacent
        # (1, 3) and (5, 7) are NOT adjacent
        intersection = max(a[0], b[0]), min(a[1], b[1])
        return intersection[0] - 1 <= intersection[1]

    # TODO: check if this works for minus strand
    intervals = sorted(set(intervals))
    combined = intervals.pop()
    multiregion = []
    while intervals:
        interv = intervals.pop()
        if adjacent(combined, interv):
            combined = union(combined, interv)
        else:  # cannot expand interval any further
            multiregion.append(combined)
            combined = interv
    multiregion.append(combined)
    return tuple(reversed(multiregion))


def isoform_display_name(s):
       """Convert clone accession ID to display friendly format"""
       return s.split('|')[0] + '-' + s.split('|')[1].split('/')[0]


def render_iso_image(orfs_to_plot, ax=None, mode='all', dname='output_isoimages',
                     intron_spacing=30, comp=300, spacing=0.5,
                     height=0.1, sp=[-1.5], show_abs_frm=False, rel_frm=False,
                     dom=False, subtle_threshold=20, mat_squares='colony'):
    """Render iso-images.
          Input:
            orfs_to_plot - list of orf_objs to plot., or set (will convert) or Isoform_Matrix
            mode - 'all'=plot skinny UTRs and thick cds exons, 'cds'=plot only thick cds exons
            blocks - plot 'exon' or 'cds' blocks
            dname - name of dir for images
            intron_spacing - padding between introns
            compression - compression (squeeze) factor
            spacing - spacing between isoform tracks
            height - height of exon blocks
            sp - list that shows the x-axis position of isoname labels (and more, depending on analysis)
            show_abs_frm - boolean, option to display absolute frm by alpha
            subtle_threshold - whether to detect and mark indiscernable subtle splicing by plotting the number corr. to num. nt difference
            mat_squares - 'colony' to plot colony images, 'color' to plot colored squares
          Output:
            pyx-generated pdf of images
        Note: assume that all orfs to plot from same genetic locus
    """

    # initialize matplotlib fig and axes objects and parameters
    max_x, min_y = 0, 0 # track abs figsize (to correctly proportion plot)
    if ax is None:
        ax = plt.gca()
    line = 0 # initial track position (y-axis)

    # process orfs to get ready for plotting
    gen_obj = grab_gen_objs_from_orfs(orfs_to_plot) # all gen_objs into a set
    verify_only_one_gen_obj_represented(gen_obj)
    verify_all_orfs_of_same_strand(orfs_to_plot)
    compress_introns_and_set_relative_orf_exon_and_cds_coords(gen_obj,
                                                orfs_to_plot, intron_spacing)
    find_and_set_subtle_splicing_status(orfs_to_plot, subtle_threshold)

    # plot gene info and header lines
    repr_orf = get_repr_orf(orfs_to_plot)
    line -= spacing/2.0
    header_line = line # save y-axis position of header, to go back and write partners
    line -= spacing/2.0

    # initialize dict mapping each orf to list of plotted rectangles
    rect_dict = {orf.name: [] for orf in orfs_to_plot}

    # workhorse loop to plot orfs and features
    # for orf in sort_orfs_by_plot_order(orfs_to_plot):
    for orf in orfs_to_plot:
        orf_name = retrieve_orf_name(orf)
        orfid = retrieve_orfid_if_exists(orf)
        label = orf_name
        # label = isoform_display_name(orf_name)
        ax.text(sp[0], line, label, va='center') # write-out orfname
        intron_start, intron_end, intron_line = get_intron_plot_coordinates(orf, comp, line, height)
        ax.add_line(mlines.Line2D([intron_start, intron_end], [intron_line, intron_line], lw=1.5, color='k', ls='--', dashes=(1,1.2), zorder=1))
        max_x = update_figure_range(max_x, intron_end)

        # render thin exon blocks
        for exon in orf.exons:
            x, y, blength, bheight = get_exon_plot_coordinates(exon, comp, height, line)
            col = get_orf_color(orf)
            alpha_val = get_abs_frm_specific_alpha_values_if_option_set(exon, show_abs_frm)
            if mode == 'all':
                ax.add_patch(mpatches.Rectangle([x, y], blength, bheight, lw=1, ec='k', fc='w', zorder=1.5, joinstyle='round')) # base layer so 'alpha' diff for 3 diff frm isn't show-through
                ax.add_patch(mpatches.Rectangle([x, y], blength, bheight, lw=1, ec='k', fc=col, zorder=2, joinstyle='round', alpha=0.2))
            # add subtle splice (delta) amounts, if option turned on
            # first, make sure the exon contains a (coding) cds object
            if exon.cds:
                delta_start, delta_end = retrieve_subtle_splice_amounts(exon.cds)
                if delta_start:
                    ax.text(x, y+height*2, delta_start, va='bottom', ha='left', size='x-small')
                if delta_end:
                    ax.text(x+blength, y+height*2, delta_end, va='bottom', ha='right', size='x-small')
            # render cds blocks, if exists
            if exon.cds:
                x, y, blength, bheight = get_exon_plot_coordinates(exon.cds, comp, height, line, cds=True)
                ax.add_patch(mpatches.Rectangle([x, y], blength, bheight, lw=1, ec='k', fc='w', zorder=3, joinstyle='round')) # base layer so 'alpha' diff for 3 diff frm isn't show-through
                exon_rect = mpatches.Rectangle([x, y], blength, bheight, lw=1, ec='k', fc=col, zorder=4, joinstyle='round', alpha=alpha_val)
                ax.add_patch(exon_rect)
                rect_dict[orf.name].append(exon_rect)
            max_x = update_figure_range(max_x, (x+blength))
            # render features (domains or isrs)
            # TODO: - debug, now that I redesigned features, exon doesn't have 'maps'
            # for m in exon.maps:
            #     fx, flen = get_feature_ranges(m, comp)
            #     if m.feat.cat == 'isr': # isrs as hatch marks
            #         ax.add_patch(mpatches.Rectangle([fx, y], flen, bheight, lw=0, fill=None, hatch='////', zorder=4))
            #     elif dom: # dom as orange blocks
            #         ax.add_patch(mpatches.Rectangle([fx, y], flen, bheight, lw=0, fc='darkorange', zorder=3, alpha=0.7))
            #         plt.text(fx, y+height*2, m.feat.name, va='bottom', size=6)
        line -= spacing
    line += spacing

    # ax.axis('off') FIXME: uncomment
    ax.axis('image')
    max_y = line - spacing
    max_x = max_x + abs(sp[0]) # adjust for sp[0] plotted

    return rect_dict


def render_pair_align_image(align_group):
    """Render iso-image for a pairwise alignment of ORFs.
    """

    # TODO: instead of modifying exon rects, draw new rects over segments of exons w/ rel frameshift
    rect_dict = render_iso_image([align_group.anchor_orf, align_group.other_orf])
    rects = rect_dict[align_group.other_orf.name]
    for exon, rect in zip(align_group.other_orf.exons, rects):
        rect.set_hatch("xx")

def grab_gen_objs_from_orfs(orfs: Iterable['ORF']) -> Set['Gene']:
    return set([orf.gene for orf in orfs])


def verify_only_one_gen_obj_represented(gen_objs: Set['Gene']):
    """Given a group of orfs to plot, verify only one gen_obj represented."""
    if len(gen_objs) > 1:
        names = [gen_obj.name for gen_obj in gen_objs]
        raise ValueError('Found multi gen_objs: {}'.format(','.join(names)))


def verify_all_orfs_of_same_strand(orfs: Iterable['ORF']):
    strands = set([orf.strand for orf in orfs])
    if len(strands) > 1:
        random_orf_name = list(orfs)[0].gene.name
        raise ValueError('found pos and neg orf from gene:' + random_orf_name)


def compress_introns_and_set_relative_orf_exon_and_cds_coords(gen_objs,
                                                orfs_to_plot, intron_spacing):
    """Iteratively traverse the range of the gene, squeezing intronic regions.
        In the process, all exon.rel_start and exon.rel_end are set.
        Overview:
            get lowest coord
            grab exons, (for neg. strand exons, reverse coords)
            iteratively traverse length of gene, when intron encountered it will jump to next exon
            at termination of function, all exons expected to hold .rel_start and rel_end values reflecting the compressed intron coords
    """
    first_pos = get_lowest_coord_of_genes(gen_objs)
    cur_idx = first_pos # idx that will travel across ORF/Exons to determine exonic/intronic regions and where to squeeze introns
    rel_idx = 1 # idx for each spacial pos. on final image (e.g. 1-5 for exon, 6-36 for 30 intron spacing, 37-40 for 2nd exon, etc.)
    exons = conso_exons_and_cdss_across_orfs_to_list(orfs_to_plot)
    adjust_exon_coords_to_make_ascending(gen_objs) # for neg. strand, flip coords, assign to exon.start_adj and exon.end_adj
    # compress exons (squeeze introns) by assigning new relative idx
    in_gene_body = True
    while(in_gene_body):
        exonic = False
        for exon in exons:
            if cur_idx == exon.start_adj:
                exon.rel_start = rel_idx
                exonic = True
            if cur_idx == exon.end_adj:
                exon.rel_end = rel_idx
                exonic = True
            if exon.start_adj <= cur_idx <= exon.end_adj:
                exonic = True
        if exonic:
            cur_idx += 1
            rel_idx += 1
        else:
            # find the next absolute index greater than cur_idx
            next_idx = 100000000000000 # impossibly big idx
            for exon in exons:
                if cur_idx < exon.start_adj < next_idx:
                    next_idx = exon.start_adj
            if next_idx == 100000000000000:
                # end of the orfs, break out of loop
                in_gene_body = False
            else:
                # set cur_idx to next lowest position (closest exon start position)
                cur_idx = next_idx
                rel_idx += intron_spacing
    update_orf_start_end_relative_coords_from_exon_coords(orfs_to_plot)

def get_lowest_coord_of_genes(gen_objs):
    """Get the left-most position coordinate.
    If neg. strand, then all coords reversed and first position is 0 or 1.
    Input:
    gen_objs - set of gen_objs"""
    repr_gen_obj = next(iter(gen_objs))
    #TODO: - deal with cases where overlapping gene on diff. strands (rare though)
    if repr_gen_obj.strand == '+':
        # don't trust the gen_obj start b/c sometimes combine gc and pb gene and they have different starts
        lowest_start = 100000000000000
        #TODO: - ensure it is ok that we only look at exon coords (not cds) coords to derive lowest coord
        exons = gather_exons_and_cds_of_genes(gen_objs)
        for exon in exons:
            if exon.start < lowest_start:
                lowest_start = exon.start
        first_pos = lowest_start
    elif repr_gen_obj.strand == '-':
        first_pos = 0 # all exon start/end coords reversed (from 100-150, 200-300 to 1-101, 151-201)
    return first_pos

def gather_exons_and_cds_of_genes(gen_objs):
    """But all exons and cds across gen_objs into a list."""
    all_exons = set()
    for gen_obj in gen_objs:
        for orf in gen_obj.orfs:
            for exon in orf.exons:
                all_exons.add(exon)
                if exon.cds:
                    all_exons.add(exon.cds)
    return all_exons

def conso_exons_and_cdss_across_orfs_to_list(orfs_to_plot, cds_only=False):
    """Put all exons belonging to a group of ORFs into one list.  Return list."""
    exons = []
    for orf in orfs_to_plot:
        if not cds_only:
            exons.extend(orf.exons)
        exons.extend(orf.cdss)
    return exons

def adjust_exon_coords_to_make_ascending(gen_objs):
    """For negative and positive strand, create 'adjusted' start/end coords of each exon for plotting purposes.
    If negative strand, convert exon absolute start/end coords to inverted coords.
    e.g. 2-5, 10-13, 20-25 will become 1-6, 13-16, 21-24
    Adjusted values saved to start_adj and end_adj"""
    #TODO: - should I add method to in this func to check all genes are same strand?, as of now it is in another place
    exons = gather_exons_and_cds_of_genes(gen_objs)
    g = list(gen_objs)[0] # repr gen_obj in set
    if g.strand == '-':
        top_coord = get_highest_coord_among_exons(exons)
        # convert exon start/end to relative coords
        adj_amt = top_coord + 1
        for exon in exons:
            exon.start_adj = abs(exon.end - adj_amt) # reverse coord
            exon.end_adj = abs(exon.start - adj_amt)
    elif g.strand == '+':
        for exon in exons:
            exon.start_adj = exon.start
            exon.end_adj = exon.end

def update_orf_start_end_relative_coords_from_exon_coords(orfs_to_plot):
    # update orf_obj relative start/end, based on newly assigned exon start/ends
    for orf in orfs_to_plot:
        orf.set_rel_start_and_end()


#############################################################################


#WISH - function that determines whether size of jutted out portion of a subtle splice
# is < threshold amount (size visually on the viewable pdf), and if it's smaller than threshold
# to then plot the delta # nt text
def find_and_set_subtle_splicing_status(orfs_to_plot, threshold=15):
    # !!! 191015 changed the subtle splicing status to be assoc. with exon (b/c plotting 6K w new isoacc)
    """Traverse all start/end coords for cdss and set subtle splice status.
    Subtle splice status stored under pop-up attribute 'cds.*_has_subtle_splice'
    e.g.
        cds.start_has_subtle_splice = 12
        cds.end_has_subtle_splice = 0
    If cds.*_has_subtle_splice is 0, then there is no subtle splice.
    Note - if negative strand, start/end are reversed
    """
    cdss = conso_exons_and_cdss_across_orfs_to_list(orfs_to_plot, cds_only=True)
    coords = cds_ranges_to_set(cdss)
    for cds in cdss:
        # determine if there is a delta between start and start of other cds
        # if so, set pop-up attr 'start_subtle_splice'
        st_delta = 0 # init. to nothing, set to smallest delta when comparing cds start to all other coords in orf group
        en_delta = 0
        for coord in coords:
            st, en = coord
            # first, assign delta for cds start position
            if st <= cds.start <= en:
                if cds.start != st:
                    delta = cds.start - st # num nt in subtle splice gap
                    if not st_delta:
                        st_delta = delta
                    else:
                        if delta < st_delta:
                            st_delta = delta
            # next, assign for cds end position
            if st <= cds.end <= en:
                if cds.end != en:
                    delta = en - cds.end # num nt in subtle splice gap
                    if not en_delta:
                        en_delta = delta
                    else:
                        if delta < en_delta:
                            en_delta = delta
        strand = orfs_to_plot[0].strand # grab strand, if neg. strand, then need to reverse st/en coords
        if strand == '-':
            st_delta, en_delta = en_delta, st_delta
        if st_delta > threshold:
            cds.start_subtle_splice = None
        else:
            cds.start_subtle_splice = st_delta
        if en_delta > threshold:
            cds.end_subtle_splice = None
        else:
            cds.end_subtle_splice = en_delta

def get_repr_orf(orfs_to_plot):
    repr_orf = orfs_to_plot[0] # representative orf (to extract gene info)
    return repr_orf


def retrieve_orf_name(orf: 'ORF') -> str:
    # retrieve printable orf_name
    # 6KOC-specific trim of orf.name (to exclude unnecessary contig_acc), e.g. MAX_p4_g05.TRINITY_etc -> MAX_p4_g05
    orf_name = orf.name
    if orf.cat == 'CL':
        if orf.src == '6KOC':
            orf_name = orf.name.split('.')[0]
    return orf_name

def retrieve_orfid_if_exists(orf: 'ORF') -> str:
    """If exists, retrieve orfid and return as a string.  orfid only exists for
    CL_ORF objects, not GC_ORF objects.  But sometimes GC_ORF objects are input
    to make GC+CL-containing tracks.
    """
    if hasattr(orf, 'orfid'):
        return orf.orfid
    else:
        return ''

def get_intron_plot_coordinates(orf, comp, line, height):
    # derive the intron start and end
    verify_orf_has_rel_start_end(orf)
    intron_start = (orf.rel_start)/comp
    intron_end = (orf.rel_end - orf.rel_start)/comp + intron_start
    intron_line = line + height/2.0
    return intron_start, intron_end, intron_line


def update_figure_range(max_x, intron_end):
    # for matplotlib xlim, track farthest right position
    if max_x < intron_end:
        max_x = intron_end
    return max_x

def set_size_of_ax(w,h, ax=None):
    """ w, h: width, height in inches """
    if not ax: ax=plt.gca()
    l = ax.figure.subplotpars.left
    r = ax.figure.subplotpars.right
    t = ax.figure.subplotpars.top
    b = ax.figure.subplotpars.bottom
    figw = float(w)/(r-l)
    figh = float(h)/(t-b)
    ax.figure.set_size_inches(figw, figh)


def get_exon_plot_coordinates(block, comp, height, line, cds=False):
    # derive the exon coordinates
    # block - exon or cds
    # cds - plotting cds
    # set thickness of blocks based on if it is a CDS or exon
    if not cds:
        n = 1
        block_length = block.len
    else:
        n = 2
        block_length = len(block.exon.aa_seq)
    x = block.rel_start/comp
    y = line
    if cds:
        y -= height/2.0
    blength = block_length/comp # block length
    bheight = height * n # block height
    return x, y, blength, bheight


def get_orf_color(orf):
    # define orf color
    col_cats = {'GC':'steelblue', 'CL':'mediumseagreen', 'PB':'indianred',
                '4K':'mediumpurple', '6K':'mediumpurple'}
    if hasattr(orf, "src"):
        return col_cats[orf.src]
    else:
        return col_cats[orf.cat]


#TODO: - make alpha values more extreme
def get_abs_frm_specific_alpha_values_if_option_set(exon, render_abs_frm):
    # if abs_frm option on, set alpha values at three shades based on abs. frm
    if render_abs_frm:
        afrm = exon.abs_frm
        if afrm == 2:
            alpha_val = 1
        elif afrm == 1:
            alpha_val = 0.45
        elif afrm == 0:
            alpha_val = 0.15
        else:
            raise UserWarning('unknown afrm')
    else:
        alpha_val = 1
    return alpha_val


def retrieve_subtle_splice_amounts(exon_cds):
    """Return text label for subtle splice gaps (deltas)."""
    delta_label_start = None
    delta_label_end = None
    if exon_cds.start_subtle_splice:
        delta_label_start = str(exon_cds.start_subtle_splice) + 'nt'
    if exon_cds.end_subtle_splice:
        delta_label_end = str(exon_cds.end_subtle_splice) + 'nt'
    return delta_label_start, delta_label_end


def get_feature_ranges(m, comp):
    """Relative to the exon, get the domain ranges."""
    fx = x - 1/comp + m.exon_start/comp # feature x-pos
    flen = (m.exon_end - m.exon_start)/comp # feature length
    return fx, flen





def dist(x1, y1, x2, y2):
    # calc. distance b/t two points
    # prevent clashes b/t mut. lollipops
    return ((y2-y1)**2 + (x2-x1)**2)**0.5


def find_nonoverlapping_height(mx, my, plotted, pad=0.1):
    # x - horizontal position of mut
    # y - (candidate) vertical position of mut
    # plotted - list of x,y coordinates already plotted
    # pad - minimum tolerated distance between lollipops
    # mx - mutational x-position
    # my - mutational y position
    no_overlap_found = True
    mut_space = 0
    while no_overlap_found:
        for x, y in plotted:
            if dist(mx, my, x, y) < pad:
                my += pad
                mut_space += pad
                break
        else:
            no_overlap_found = False
    return mut_space


def get_disease_acronym(dis_descr):
    # dis_descr - description line of disease
    first_chars = [x[0] for x in dis_descr.split()]
    acro = ''.join(first_chars).upper()
    return acro



def get_highest_coord_among_exons(exons):
    """Given a list or set of exon objs, return highest coord value."""
    top_coord = 0
    for exon in exons:
        if exon.end > top_coord:
            top_coord = exon.end
    return top_coord












def cds_ranges_to_set(list_of_cds_objs):
    """Given a list of cds_objs, put start/end coord to list.
    e.g. set([100,150], [200, 250])"""
    coords = set()
    for cds in list_of_cds_objs:
        coord = (cds.start, cds.end)
        coords.add(coord)
    return coords




def sort_orfs_by_plot_order(orfs_to_plot):
    # plot order by 'ordinal' (GC first, then CL), then by isoname alpha
    orfs_to_plot = sorted(orfs_to_plot, key=operator.attrgetter('plot_ord', 'name'))
    return orfs_to_plot




def verify_orf_has_rel_start_end(orf):
    # verify that the orf was updated with rel_start and rel_end
    if orf.rel_start == None or orf.rel_end == None:
        raise UserWarning('orf.rel_start/end = None:', orf.name, orf.get_exon_ranges()) # found issue when both pos and neg strand orf part of gene


#TODO: - check for accuracy of the function
def get_mutation_coords(pos, comp, mut_adjust, line, plotted):
    """Determine mutation coordinates.  Return coord. for plotting of lollipops.
    """
    mstart = x - 1/comp + pos.ridx/comp
    mline = line + mut_adjust
    mut_space = 0
    mut_space = find_nonoverlapping_height(mstart, mline, plotted)
    plotted.append((mstart, mline + mut_space)) # track already-plotted coordinates to avoid clashing for subsequent lollipops
    return mstart, mline, mut_space


def get_mut_color(mut):
    """Determine disease mutation color."""
    if mut.cat == 'Patho':
        col = 'r'
    elif mut.cat == 'Benign':
        col = 'g'
    else:
        raise UserWarning('valid mut.cat not found:' + mut.dis)
    return col


def update_disease_acronym_dictionary(mut, dis_acro_dict):
    """Add disease info to dict. for disease legend upon final write-out.
    In the process, 'dis_acro_dict' gets updated in-place.
    """
    dis_acro = get_disease_acronym(mut.dis) # print(out acronym of disease)
    dis_acro_dict[dis_acro] = mut.dis # add acronym key for legend

def is_orfs_to_plot_a_list_or_set(orfs_to_plot):
    """Determine if the passed in argument 'orfs_to_plot' is a set or list.
    If not, then it is an Isoform_Matrix.
    """
    arg_type = type(orfs_to_plot)
    if arg_type == set or arg_type == list:
        return True
    return False

def extract_orfs_from_isoform_matrix(isoform_matrix):
    """Input is an isoform_matrix object.  Therefore, need to pull out orfs.
    """
    mat_obj = isoform_matrix
    orfs_to_plot = list(mat_obj.orfs)
    return mat_obj, orfs_to_plot

def does_iso_ptr_pairset_have_at_least_one_pos_ppi(ppi, db_orf):
    """Determine if the iso-ptr that makes up the ppi, that there is at least
    one positive ppi involving the gene/iso-ptr set.
    """
    ad_orf = ppi.ad
    ad_orf_ppis = ad_orf.ppis_as_ptr
    has_pos = False
    for ppi in ad_orf_ppis:
        if ppi.db.gene.name == db_orf.gene.name:
            if ppi.final_score == '1':
                has_pos = True
    return has_pos

def create_ordered_list_of_orfid_ppi_pairs(orf, mat_gene):
    """From an orf, extract the linked ppis and put into a list of orfid/ppi
    pairs.
    Note: skips over cases of no positive PPI for the partner (i.e. all-neg partner)
    Input:
        orf - representative orf used to extract set of PPIs for that iso-matrix
        mat_gene - genename of the 'anchor', gene forming basis of iso-matrix*
            *need to track of this because sometimes orf has linked ppi, but
            because part of the prey-side of another iso-matrix (e.g. AES ppi
            as DB, but also as AD for lots of other interactions)
    Output:
        ordered_ppis - list of list with orf_id, ppi_obj pairs
    """
    ordered_ppis = [] # [[orf_id, ppi_obj], [orf_id, ppi_obj]]
    ppis = orf.ppis
    for ppi in ppis:
        if does_iso_ptr_pairset_have_at_least_one_pos_ppi(ppi, orf):
            if ppi.db.gene.name == mat_gene:
                orfid = ppi.ad.orfid
                pair = [orfid, ppi]
                ordered_ppis.append(pair)
    ordered_ppis = sorted(ordered_ppis)
    return ordered_ppis

def extract_genenames_from_ordered_ppi_list(ordered_ppis):
    """From a list of [<orfid>, <ppi_obj>] pairs, return a list of genenames."""
    genenames = []
    for pair in ordered_ppis:
        orfid, ppi_obj = pair
        genename = ppi_obj.ad.gene.name
        genenames.append(genename)
    return genenames

def extract_orfnames_from_ordered_ppi_list(ordered_ppis):
    """From a list of [<orfid>, <ppi_obj>] pairs, return a list of genenames."""
    orfnames = []
    for pair in ordered_ppis:
        orfid, ppi_obj = pair
        orfname = ppi_obj.ad.name
        orfnames.append(orfname)
    return orfnames

def extract_genenames_and_orfnames_from_ordered_ppi_list(ordered_ppis):
    name_pairs = []
    for pair in ordered_ppis:
        orfid, ppi_obj = pair
        orfname = ppi_obj.ad.name
        genename = ppi_obj.ad.gene.name
        name_pair = [genename, orfname]
        name_pairs.append(name_pair)
    return name_pairs

def generate_image_array_and_plot_coords(img_path, right_tmp, line_tmp, cpad):
    """Given the path of the colony image and top-left coordinates to plot
    the grid of 9 colony images, return the image_array (needed by matplotlib
    for plotting) and computed coordinates.
    """
    im_array = plt.imread(img_path)
    extent = (right_tmp, right_tmp+cpad, line_tmp, line_tmp+cpad) # left, right, bottom, top
    return im_array, extent

def determine_ppi_result_square_color(ppi_obj):
    """Given option to plot colored squares (black, white, etc.), determine
    which color square to plot."""
    color_map = {0:'white', 1:'black'}
    fn_score = ppi_obj.final_score
    if fn_score.isdigit():
        fn_score_as_int = int(fn_score)
    else:
        fn_score_as_int = -1
    if fn_score_as_int in color_map:
        color = color_map[fn_score_as_int]
        return color
    else:
        return None

def print_iso_char_image(orf1, orf2, scale=30):
    # helper functions
    def is_pos_in_orf(orf, i):
        for exon in orf.exons:
            for pos in exon.chain:
                if i == pos.coord:
                    return True
        return False

    def return_cat(in_orf1, in_orf2):
        if in_orf1 and in_orf2:
            return 'C'
        elif in_orf1:
            return '1'
        elif in_orf2:
            return '2'
        else:
            return '0'

    def jump_index(orf1, orf2, i):
        # set index to left-most-closest coord. position
        smallest_coord = 1000000000000
        for exon in orf1.exons + orf2.exons:
            for pos in exon.chain:
                if pos.coord > i and pos.coord < smallest_coord:
                    smallest_coord = pos.coord
        return smallest_coord

    def split_blocks(in_string):
        # split basd on contiguous char.
        split_string = [''.join(v) for k, v in groupby(in_string)]
        return split_string

    def downsample_block(size, scale):
        # downsample size of exon block by amount
        small_size = int(math.ceil(float(size)/scale))
        return small_size

    def convert_blocks_to_iso_chains(blocks, scale):
        first = ''
        second = ''
        for block in shrunk_blocks:
            for b in block:
                if b == 'C':
                    first += 'X'
                    second += 'X'
                elif b == '1':
                    first += 'X'
                    second += '-'
                elif b == '2':
                    first += '-'
                    second += 'X'
                elif b == '0':
                    first += '-'
                    second += '-'
        return first, second

    # print(out an intron-squeezed isoform representation)
    # start with lowest coord, ascend and set category
    i = 0 # current coordinate
    i_orf1 = 0 # current index of exon obj for orf1
    i_orf2 = 0 # current index of exon obj for orf2
    # first, find the lowest abs. coordinate between two orfs
    coord1 = orf1.exons[0].chain[0].coord
    coord2 = orf2.exons[0].chain[0].coord
    if coord1 == coord2:
        i = coord1
        cat = 'C'
    elif coord1 < coord2:
        i = coord1
        cat = '1'
    elif coord2 < coord1:
        i = coord2
        cat = '2'
    chain = ''
    chain += cat
    for j in range(1000):
        i += 1
        in_orf1 = is_pos_in_orf(orf1, i)
        in_orf2 = is_pos_in_orf(orf2, i)
        cat = return_cat(in_orf1, in_orf2)
        chain += cat
        if cat == '0':
            i = jump_index(orf1, orf2, i) - 1
        if i == 1000000000001:
            break
    blocks = split_blocks(chain)
    shrunk_blocks = []
    for block in blocks:
        cat = block[0]
        size = downsample_block(len(block), scale)
        shrunk_blocks.append(cat*size)
    first, second = convert_blocks_to_iso_chains(blocks, scale)
    print(first)
    print(second)
