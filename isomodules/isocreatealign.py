#!/usr/bin/env python
# title           :isocreate.py
# description     :Functions to create web of isoform objects.
# author          :Gloria Sheynkman
# date            :May 2nd, 2019
# version         :1
# python_version  :2.6.6
# ==============================================================================

from . import isoclass
from collections import Counter
from . import isoalign
from . import isogroup
import itertools
from itertools import groupby
import re


# *****************************************************************************
# create and map alignemnts between orfs, based on genome align
def create_and_map_splice_based_align_obj(orf_pairs, abacus=True):
    """Based on genome-centric align, create aln_obj for res/cds/orf.
       In the process, 'call' categories of alignments.
       Returns a list of grps, each grp has orf pair and linked aln_objs.
       Gene_dict is modified in-place.

       Note - In the process of creating alignment object, sets the relative
              frame for position and residue objects. This is needed to determine
              the residue match category.

       Input:
        orf_pairs - list of orf pairs (assume first orf is 'ref')
        abacus - abacus the nt overhangs at junctions

       Example of why abacusing needed:
        1----23
        1-23123
       Algorithm takes first nt of codon and grabs other-res with overlap of
       2 or 3 nt. So in this e.g., will take the second res of alt. isoform.
       Therefore, do abacusing:
        ----123
        -123123
    """
    grps = []  # group of orf pairs with populated aln_objs
    for pair in orf_pairs:
        orf1, orf2 = pair  # orf1 = anchor_orf, orf2 = other_orf
        set_coords(pair, abacus)
        all_coords = get_list_of_coords(pair)
        orf1_coords = make_coord_dict(orf1)
        orf2_coords = make_coord_dict(orf2)
        # create orf-pair-flushed chain of res_obj {orf1.name:[], orf2.name:[]}
        # for res_obj that don't exist, create empty res_obj to fill in
        chains = make_orf_pair_aligned_aa_chain(all_coords, pair, orf1_coords,
                                                orf2_coords)
        set_rfrm_of_pos_in_orf(all_coords, orf1_coords, orf2_coords, pair)
        set_rfrm_of_res_in_chain(chains)

        # make a grp_obj of the two orfs
        grp = isogroup.PairwiseAlignmentGroup(orf1, orf2)
        grps.append(grp)
        # instantiate the 'full' alignment (between orfs)
        alnf = isoalign.AlignmentFull(grp, orf1, orf2)
        grp.alnf = alnf

        # make and link 'residue' aln_objs (via alnr obj)
        for res1, res2 in zip(chains[orf1.name], chains[orf2.name]):
            match = get_the_match_type_of_the_two_residues(res1, res2)
            # create a residue-level alignment
            aln_obj = isoalign.AlignmentResidue(alnf, match, res1, res2)
            alnf.chain.append(aln_obj)

        # make and link 'block' alnb_objs
        ranges = get_ranges_of_contiguous_blocks_w_same_match_type(alnf.chain)
        for match_cat, start, end in ranges:
            alnr_chain = alnf.chain[start: end]
            alnb = isoalign.AlignmentBlock(match_cat, alnf, alnr_chain)

        # find 'protein-centric' blocks, make and link 'pblock' alnpb_obj
        # alnf.blocks represent splice-based effects
        #  For example: IDMDMFID
        # alnf.pblocks represent protein-based effects
        #  For example: SMDMS (drived from e.g. above)
        # ---
        # first, run through and pinpoint cases where subst. is actually
        # a match (e.g. MA--V, --MAV)
        # output new blocks e.g. IDMDM -> MMMDM -> ['MMM', 'D', 'M']
        # note - mismatch at exon edge included as part of ins (I) and del (D)
        full_block_string, split_blocks = merge_blocks_based_on_protein_effect(alnf)
        i = 0
        new_full_block_string = ''
        for block_string in split_blocks:
            cat = get_the_category_of_the_block(block_string)
            # grab correpsonding alnb_objs
            start, end = i, i + len(block_string)
            i = end
            # TODO: - 200607 had nonASCII character, check code below
            alnbs = alnf.blocks[start: end]
            # temporary alnpb to find if same residue
            alnpb = isoalign.AlignmentProteinBlock(cat, alnf, alnbs)
            if cat == 'S' and alnpb.aa1 == alnpb.aa2:
                new_block_string = 'M' * len(block_string)
            else:
                new_block_string = block_string
            new_full_block_string += new_block_string
        new_split_blocks = split_blockstring_w_repeating_Ms(new_full_block_string)
        # second, find 'protein-centric' blocks, make and link 'pblock' alnpb_obj
        # i = 0
        # for block_string in new_split_blocks:
        #     cat = get_the_category_of_the_block(block_string)
        #     # grab correpsonding alnb_objs
        #     start, end = i, i + len(block_string)
        #     i = end
        #     alnbs = alnf.blocks[start: end]
        #     alnpb = isoalign.AlignmentProteinBlock(cat, alnf, alnbs)

        # define and create 'subblock' aln_objs
        # note - functions below also used for clustal-based alignment
        subblock_ranges = determine_subblock_ranges(alnf.chain)
        for start, end in subblock_ranges:
            alnsb_chain = get_subset_of_alnr_based_on_range(alnf, start, end)
            cds1 = get_cds_mapped_to_alnr_chain(alnsb_chain, 1)
            cds2 = get_cds_mapped_to_alnr_chain(alnsb_chain, 2)
            alnsb = isoalign.AlignmentSubblock(alnf, cds1, cds2, alnsb_chain)
    return grps


def set_coords(pair, elected_abacus_option):
    """Determine if there are one-nt overhangs of codons acros junctions.
    If so, "abacus" the nt to form a contiguous codon.  By set of pos_obj
    attribute coord_v2.
    """
    for orf in pair:
        for res in orf.res_chain:
            # TODO: - can we move 'coord_v2' attr into the Position class?
            res.p1.coord_v2 = res.p1.coord
            res.p2.coord_v2 = res.p2.coord
            res.p3.coord_v2 = res.p3.coord
            if elected_abacus_option:
                if orf.strand == '+':
                    if (res.p1.coord + 3) < res.p2.coord:
                        res.p1.coord_v2 = res.p2.coord - 1
                    elif (res.p2.coord + 3) < res.p3.coord:
                        res.p3.coord_v2 = res.p2.coord + 1
                elif orf.strand == '-':
                    if (res.p1.coord - res.p2.coord) != 1:
                        res.p1.coord_v2 = res.p2.coord + 1
                    elif (res.p2.coord - res.p3.coord) != 1:
                        res.p3.coord_v2 = res.p2.coord - 1

def get_list_of_coords(orf_pair):
    """Return a list of absolute coordinates of pos. obj in two orfs.
       If + strand, list is ascending.
       If - strand, list is descending.
       Note - only considering protein-coding pos_objs.
    """
    coords = set()
    for orf in orf_pair:
        for pos in orf.chain:
            if pos.res:
                coords.add(pos.coord_v2)
    if orf.strand == '+':
        return sorted(list(coords))
    elif orf.strand == '-':
        return sorted(list(coords), reverse=True)

def make_coord_dict(orf):
    """Make and return a dict. of coord->pos_obj of all pos_obj of
       orf. Assumes that pos.coord_v2 are populated.
       Note - only considering protein-coding pos_objs.
    """
    coords = {}  # coord_v2->pos_obj
    for pos in orf.chain:
        if pos.res:
            coords[pos.coord_v2] = pos
    return coords

def make_orf_pair_aligned_aa_chain(all_coords, pair, orf1_coords, orf2_coords):
    """Given the coord ranges and two orf_obj, get the genome-spaced (i.e.,
       exactly-acoordioned, i.e., flushed aa chains) res_objs.
    """
    orf1, orf2 = pair
    chains = {orf1.name: [], orf2.name: []}  # orf-pair-aligned res_obj chains
    first_coord = all_coords[0]
    last_coord = all_coords[-1]
    coord = first_coord
    while(coord != last_coord):
        anchor_res = grab_an_overlapping_residue(coord, pair)
        if anchor_res:
            anchor_orf = anchor_res.orf
            other_orf = get_other_orf(anchor_orf, pair)
            other_res = get_overlapping_residue_if_exists(anchor_res, pair, orf1_coords, orf2_coords)
            if not other_res:
                other_res = isoclass.EmptyResidue(isoclass.EmptyCDS())
            chains[anchor_orf.name].append(anchor_res)
            chains[other_orf.name].append(other_res)
            coord = get_latest_first_nt_of_codon_of_two_residues(anchor_res, other_res)
        else:
            anchor_res = isoclass.EmptyResidue(isoclass.EmptyCDS())
        try:
            coord = get_next_coord_in_all_coords(coord, all_coords)
        except:
            print(orf1.name)
            raise UserWarning('coord does not exist:' + '\n'.join(map(str, all_coords)) + '\n\n' + str(coord))
    return chains

def grab_an_overlapping_residue(abs_coord, pair):
    """Given an absolute coordinate, return at least one residue of pair.
       Only grab residue if it corresponds to first codon.
    """
    orf1, orf2 = pair
    for res in orf1.res_chain:
        if res.p1.coord_v2 == abs_coord:
            return res
    for res in orf2.res_chain:
        if res.p1.coord_v2 == abs_coord:
            return res
    return None

def get_overlapping_residue_if_exists(res_obj, pair, orf1_coords, orf2_coords):
    """Given a res_obj, grab the overlapping residue, defined as
       an overlap of 2 or more nt of the codon.
       In some cases, check non-abacus'd coords for overlap.
    """
    orf1, orf2 = pair
    other_orf, other_coords = from_anchor_res_obj_get_other_orf_info(res_obj,
                                         orf1, orf2, orf1_coords, orf2_coords)
    codon_coords = [res_obj.p1.coord_v2,
                    res_obj.p2.coord_v2,
                    res_obj.p3.coord_v2]
    # find overlapping residue
    # criteria - if 2+ coords match pos of same residue in other orf
    # first check using abacused coords (.coord_v2), then if only
    # partial matches (1/3), check based on non-abacus'd coords
    other_res = get_overlapping_res_info(codon_coords, other_coords)
    if not other_res:
        other_res = try_grabbing_res_based_on_non_abacused_coords(res_obj, other_orf)
    return other_res

def from_anchor_res_obj_get_other_orf_info(res_obj, orf1, orf2, orf1_coords,
                                                                  orf2_coords):
    if res_obj.orf == orf1:
        other_orf = orf2
        other_coords = orf2_coords
    elif res_obj.orf == orf2:
        other_orf = orf1
        other_coords = orf1_coords
    else:
        raise UserWarning('did not localize orf_obj for res_obj')
    return other_orf, other_coords

def get_overlapping_res_info(codon_coords, other_coords):
    """Determine if the first res overlaps with res in other orf by at least
       two position objects.
    """
    other_res_objs = []
    for coord in codon_coords:
        if coord in other_coords:
            other_pos_obj = other_coords[coord]
            other_res_objs.append(other_pos_obj.res)
    cts = Counter(other_res_objs)
    overlapping_res_obj = None
    for val, ct in cts.items():
        if ct == 3 or ct == 2:
            overlapping_res_obj = val
            break
    return overlapping_res_obj


def try_grabbing_res_based_on_non_abacused_coords(res_obj, other_orf):
    """In this case, there was no overlapping res using abacus'd coords.
       It's possible that frame-shifted, exon overhang abacus'ing caused
       non-overlap.
       ex  from ZNF879-201, ZNF879-202
           [G, G, G] 54-G [179027599, 179028032, 179028033]
           [G, G, G] 65-G [179027598, 179027599, 179028032]
           *99 was abacused to downstream. *32 in second was abacused upstream.
       This caused nonoverlap.
    """
    other_res_objs = []
    other_res_obj = None
    codon_coords = [res_obj.p1.coord,
                    res_obj.p2.coord,
                    res_obj.p3.coord]
    for pos in other_orf.chain:
        if pos.coord in codon_coords:
            if pos.res:
                other_res_objs.append(pos.res)
    cts = Counter(other_res_objs)
    for val, ct in cts.items():
        if ct == 3 or ct == 2:
            other_res_obj = val
    return other_res_obj


def get_other_orf(anchor_orf, pair):
    """Given a res_obj, return the orf that is the 'other' orf not corr. to the
       orf of this res_obj.
    """
    orf1, orf2 = pair
    if anchor_orf == orf1:
        other_orf = orf2
    elif anchor_orf == orf2:
        other_orf = orf1
    else:
        raise UserWarning('res_obj not matching either orf')
    return other_orf


def get_latest_first_nt_of_codon_of_two_residues(res1, res2):
    """Given two res_obj (either of which can be None), determine the next
       absolute coord which represents the most downstream
       first-nt-of-codon of the two residues.
    """
    # initialize starting coords, strand-aware
    if res1.orf.strand == '+':
        coord1, coord2 = 0, 0
    else:
        coord1, coord2 = 1e30, 1e30
    if res_is_not_empty(res1):
        coord1 = res1.p1.coord_v2
    if res_is_not_empty(res2):
        coord2 = res2.p1.coord_v2
    if res1.orf.strand == '+':
        next_coord = max(coord1, coord2)
    elif res1.orf.strand == '-':
        next_coord = min(coord1, coord2)
    return next_coord

def res_is_not_empty(res_obj):
    if res_obj.name != '-':
        return True
    else:
        return False

def get_next_coord_in_all_coords(coord, all_coords):
    """Get the next coord in the set of nts of two orfs."""
    i = all_coords.index(coord)
    next_coord = all_coords[i + 1]
    return next_coord


def set_rfrm_of_pos_in_orf(all_coords, orf1_coords, orf2_coords, pair):
    """Set attribute pos_obj.rfrm (relative frame) for all pos_obj of two
       orfs.  Define the primary orf as the one with the longest nt chain.
       Assume primary is 0-index of pair.

       pos.rfrm set to '*' for seq. in the other_orf but not anchor_orf
    """
    orf1, orf2 = pair
    anchor_coords = orf1_coords
    other_coords = orf2_coords
    for coord in all_coords:
        if coord in anchor_coords:
            anchor_pos = anchor_coords[coord]
            anchor_pos.rfrm = 1
            if coord in other_coords:
                other_pos = other_coords[coord]
                other_pos.rfrm = calc_rfrm(anchor_pos.frm, other_pos.frm)
        else:
            other_pos = other_coords[coord]
            other_pos.rfrm = '*'

def calc_rfrm(anchor_frm, other_frm):
    if anchor_frm == other_frm and anchor_frm != None:
        return 1
    elif (anchor_frm, other_frm) in ((0, 2), (1, 0), (2, 1)):
        return 2
    elif (anchor_frm, other_frm) in ((0, 1), (1, 2), (2, 0)):
        return 3
    else:
        raise UserWarning('incompatible anchor/other frame:{},{}'.format(anchor_frm, other_frm))


def set_rfrm_of_res_in_chain(chains):
    """Set temporary rfrm attribute for residues in the aa_chain."""
    for orfname, aa_chain in chains.items():
        for res in aa_chain:
            if res_is_not_empty(res):
                # if rfrm is *,   then try to find non-* rfrm
                rfrm = str(res.p1.rfrm)
                if rfrm == '*':
                    rfrm = str(res.p2.rfrm)
                if rfrm == '*':
                    rfrm = str(res.p3.rfrm)
                res.rfrm = rfrm
            else:
                # no res_obj exists, because res in one orf but not the other
                pass


def get_the_match_type_of_the_two_residues(res1, res2):
    """Compare the AA seq. of the two residues. Return a enum descr. the
       comparison.
       Key: M=same, I=insertion, D=deletion, F=frameshift
    """
    cat = ''
    if res1.aa == res2.aa and res1.rfrm == res2.rfrm:
        cat = 'M'
    elif res1.aa == '-' and res2.aa != '-':
        cat = 'I'
    elif res1.aa != '-' and res2.aa == '-':
        cat = 'D'
    elif res1.rfrm != res2.rfrm:
        cat = 'F'
    else:
        cat = 'X'
    return cat


def get_ranges_of_contiguous_blocks_w_same_match_type(alnr_chain):
    """From a chain of aln residues, group same-match-category residues.
       Return info. on length of contigous blocks.
       Return range_info - e.g. MMMIIM -> [('M', 0, 3), ('I', 3, 5), ('M', 5, 6)]
    """
    cats = [x.cat for x in alnr_chain]
    block_info = [(k, sum(1 for i in g)) for k, g in groupby(cats)]
    # make ranges
    range_info = []  # [('M', 0, 3), ('I', 3, 5), ('M', 5, 6)]
    i = 0
    for cat, blen in block_info:
        start = i
        end = i + blen
        i = end
        range_info.append([cat, start, end])
    return range_info


def merge_blocks_based_on_protein_effect(alnf):
    """alnf.blocks represent splice-based effects. Merge blocks based on the
       protein-based effects.
       For example: IDMDMFID -> SMDMS
    """
    block_string = ''.join([alnb.cat for alnb in alnf.blocks])  # e.g. IDMDMFID
    split_blocks = filter(None, re.split('(M)', block_string))
    return block_string, split_blocks


def get_the_category_of_the_block(block_string):
    """Categorize the splice-based blocks into protein-based blocks.
       Cat changes:
        M -> M, I -> I, D -> D, any combo of I/D/F -> S (subst)
    """
    if block_string in ['M', 'I', 'D', 'DX', 'XD', 'IX', 'XI']:
        return block_string.replace('X', '')
    elif set(list(block_string)) == set(['M']):
        return 'M'
    else:
        return 'S'


def split_blockstring_w_repeating_Ms(new_full_block_string):
    """Redefined the block categories to call I/D/F (subst) tracks with same
       sequence as M (match).
       e.g. IDMDM -> MMMDM -> ['MMM', 'D', 'M']
    """
    split_blocks = filter(None, re.split('(M+)', new_full_block_string))
    return split_blocks


def determine_subblock_ranges(alnr_chain):
    """Determine the ranges corresponding to 'subblocks', defined as
       contiguous runs of same-anchor_cds/other_cds/aln_block_cat.
       For example:
        cds1 : 1122233333
        cds2 : 1122233444
        aln  : iiiiisssss
        i = identical, s = substitution
        -> ranges are 0-1, 2-4, 5-6, 7-9, because same cds1/cds2/aln
        Note - another variant of this function for clustal-align
    """
    mode_chain = []  # e.g. [(1, 1, i), (1, 1, i), etc.]
    for alnr in alnr_chain:
        cds1 = alnr.res1.cds.ord
        cds2 = alnr.res2.cds.ord
        aln = alnr.alnb.cat  # cat of aln 'block' (alnb), e.g. 'ident'
        mode_chain.append((cds1, cds2, aln))
    subblock_ranges = get_ranges_of_same_mode_subblocks(mode_chain)
    return subblock_ranges

def get_ranges_of_same_mode_subblocks(mode_chain):
    """Given a list of 'modes', return ranges of contiguous same-mode groups.
       For example: [(1,1,i),(1,1,i),(2,2,i)] -> [(0,1), (2,2)]
    """
    blens = []
    for k, g in itertools.groupby(mode_chain):
        blens.append(len(list(g)))
    starts = [sum(blens[:i]) for i in range(len(blens) + 1)][:-1]
    ends = [sum(blens[:i]) - 1 for i in range(len(blens) + 1)][1:]
    return zip(starts, ends)


def get_subset_of_alnr_based_on_range(alnf, start, end):
    return [alnf.chain[i] for i in range(start, end + 1)]


def get_cds_mapped_to_alnr_chain(alnr_chain, rank):
    """Retrieve the cds_obj the alnr are linked to.
       Based on range computed earlier, expectation is the group of alnrs
       will map to only one cds. If it maps to more than one cds, throw error.
       Note - sometimes will be mapped to no cds (because insertion/deletion)
    """
    mapped_cds = set()
    if rank == 1:
        for alnr in alnr_chain:
            mapped_cds.add(alnr.res1.cds)
    elif rank == 2:
        for alnr in alnr_chain:
            mapped_cds.add(alnr.res2.cds)
    if there_is_a_single_empty_cds(mapped_cds):
        return list(mapped_cds)[0]
    if len(mapped_cds) > 1 or len(mapped_cds) == 0:
        # there are 2 distinct cdss, remove the emptycds
        print('error - failure of one cds mapped by alnrs ' + str(alnr))
    return list(mapped_cds)[0]

def there_is_a_single_empty_cds(mapped_cds):
    """Determine if the alnrs all mapped to empty CDSs."""
    cds_names = list(set([cds.name for cds in mapped_cds]))
    if cds_names == ['-']:
        return True
    return False


def get_alnb_mapped_to_alnr_chain(alnr_chain):
    """Get the associated aln 'block' associated with the alnr of an alnsb.
       If population of alnr/alnb/alnsb was correct, expect only one alnb.
    """
    mapped_alnb = set()
    for alnr in alnr_chain:
        mapped_alnb.add(alnr.alnb)
    if len(mapped_alnb) > 1:
        pass
        # print('error - aln subblock maps to 2+ aln block ' + str(alnr.alnf))
    if len(mapped_alnb) == 0:
        print('no mapped alnb from alnr of an alnsb ' + str(alnr.alnf))
    return list(mapped_alnb)[0]
