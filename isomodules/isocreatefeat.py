#!/usr/bin/env python
# title           :isocreate.py
# description     :Functions to create features mapped to isoforms.
# author          :Gloria Sheynkman
# date            :June 17th, 2019
# version         :1
# python_version  :2.6.6
# ==============================================================================

from collections import defaultdict
from isomodules import isofeature
import itertools
from itertools import groupby

# create and map generic feature to iso_obj
# note - similar function below with "domain" flavor
def create_and_map_feature(orf, start, end, ftype, desc=None):
    """Given an anchor orf and a range on that orf (start and end AA), create
       a hierarchical feat_obj mapping (featf, featsbs, and featrs).

       ftype can be: dom, ptm, frm
    """
    featf = isofeature.FeatureFull(orf, ftype, desc=desc)
    feat_ress = orf.res_chain[int(start) - 1 : int(end)]
    for res in feat_ress:
        featr = isofeature.FeatureResidue(featf, res)
        featf.chain.append(featr)
    featr_groupings = find_groups_of_featr_in_subblocks(featf.chain)
    for (cds_ord, cds), featrs in sorted(featr_groupings.items()):
        featsb = isofeature.FeatureSubblock(featf, cds, featrs)
    return featf


# *****************************************************************************
# create and map features to iso_obj
def create_and_map_domain(orf, domain_info):
    """Given an orf and info. on a domain mapped to its AAs, create a
       hierarchical feat_obj mapping.
       Domain objects contain the domf, domsbs, and domrs.
       Input:
        orf - orf_obj
        domain_info - [pfam_acc, pfam_name, cat (DBD, other), eval,
                       1-base_start, 1-base_end]
      Return (optional):
       featf (domain object)
    """
    # TODO - issue is that if this is run twice on same orf,
    #        multiple identically-named domains will map to orf
    pfam, name, cat, eval, start, end = domain_info
    featf = isofeature.DomainFull(orf, cat, pfam, name, eval)
    domain_ress = orf.res_chain[int(start) - 1 : int(end)]
    for res in domain_ress:
        featr = isofeature.DomainResidue(featf, res)
        featf.chain.append(featr)
    featr_groupings = find_groups_of_featr_in_subblocks(featf.chain)
    for (cds_ord, cds), featrs in sorted(featr_groupings.items()):
        featsb = isofeature.FeatureSubblock(featf, cds, featrs)
    return featf

def find_groups_of_featr_in_subblocks(featrs):
    """Given a list of featr (repr. domain mapped to orf), group the featrs
       into those mapping to the same cds.

       Return:
        {(cds_ordinal, cds_obj) -> [featr, featr, featr]}
       cds_ordinal allows for sorting of keys so featsb can be populated
       into featf in upstream-to-downstream order
    """
    sb_grps = defaultdict(list)
    for featr in featrs:
        cds_ord = featr.res.cds.ord
        cds = featr.res.cds
        acc = (cds_ord, cds)
        sb_grps[acc].append(featr)
    return sb_grps



# NOTE - frame feature type, as implemented only exists for the other orf
# consider editing code so that a featf (frmf) is made for the anchor orf too
# *****************************************************************************
# create and map frm feature to iso_obj
def create_and_map_frame_objects(grp):
    """Given a group which holds two orfs which are linked through an alignment
       object, create a frame object (as a feature) covering the full length
       of the other (not anchor/repr, because frame would be 1) orf.
       Note - assumes that the res.rfrm attribute is pre-populated from
              creation of the alnf which exactly corresponds to the same grp.
    """
    orf2 = grp.other_orf
    fstr = ''.join([res.rfrm for res in orf2.res_chain])
    frame_blocks = [''.join(g) for _, g in itertools.groupby(fstr)]
    infer_list, frm_list = infer_frame_of_non_overlapping_seq(frame_blocks)
    frmf = isofeature.FrameFull(orf2, grp)
    grp.frmf = frmf
    for i, frm in enumerate(frm_list):
        is_inferred = infer_list[i] # 0 or 1
        res = orf2.res_chain[i]
        infer_status = 'direct' if not is_inferred else 'indirect'
        frmr = isofeature.FrameResidue(frmf, res, frm, is_inferred)
        frmf.chain.append(frmr)
    ranges = get_ranges_of_contiguous_blocks_w_same_frm(frmf.chain)
    for match_cat, start, end in ranges:
        frmrs, frmb_cat = extract_the_frmb_info(frmf, start, end)
        frmb = isofeature.FrameBlock(frmf, frmrs, frmb_cat)

def infer_frame_of_non_overlapping_seq(frame_blocks):
    """Infer the frame of sequence regions of alt. iso. that are not in ref.
       Algorithm - look to the upstream-most overlapping segment. If doesn't
       exist (b/c is at N-terminus), look to downstream.
    """
    inferred_frame_blocks = []
    infer_status_chain = []
    for i, block in enumerate(frame_blocks):
        blen = len(block)
        if i == 0 and '*' in block:
            inferred_frm = frame_blocks[1][0]
            infer_status = True
        elif '*' in block:
            inferred_frm = frame_blocks[i - 1][0]
            infer_status = True
        else:
            # don't need to infer frame b/c propogated from ref.
            inferred_frm = frame_blocks[i][0]
            infer_status = False
        infer_status = infer_status * 1  # convert to no. for writeout
        inferred_frame_blocks.append(inferred_frm * blen)
        infer_status_chain.extend([infer_status] * blen)
    frm_list = list(''.join(inferred_frame_blocks))
    infer_list = list(''.join(map(str, infer_status_chain)))
    return infer_list, frm_list

def create_list_of_frms_and_ifrms(inferred_blocks):
    ifrm_str = ''.join(inferred_blocks)
    return ifrm_str

def get_ranges_of_contiguous_blocks_w_same_frm(frmr_chain):
    """From a chain of frm residues, group same-frm-category residues.
       Return info. on length of contigous blocks.
       Return range_info - e.g. 111221-> [('1', 0, 3), ('2', 3, 5), ('1', 5, 6)]
       frmr.cat can be 1, 2, or 3
    """
    cats = [x.cat for x in frmr_chain]
    block_info = [(k, sum(1 for i in g)) for k, g in groupby(cats)]
    # make ranges
    range_info = []  # [('1', 0, 3), ('2', 3, 5), ('1', 5, 6)]
    i = 0
    for cat, blen in block_info:
        start = i
        end = i + blen
        i = end
        range_info.append([cat, start, end])
    return range_info

def extract_the_frmb_info(frmf, start, end):
    """Get the info. needed to instantiate a frmb object.
       In the process, ensure that all extracted frmrs are of the same frame (cat).
       Return a chain of same-frame frmrs, the cat (frame '1', '2', or '3')
    """
    frmrs = frmf.chain[start: end]
    frmb_cat = QC_filter_only_one_cat_in_frmb(set([frmr.cat for frmr in frmrs]))
    return frmrs, frmb_cat

def QC_filter_only_one_cat_in_frmb(frm_cats):
    """Previously split the frmrs based on same-frame blocks. QC check that
       only one frame (frmr.cat) is represented.
       Input frm_cats is a set of frm_cats
    """
    if len(frm_cats) != 1:
        raise Warning('Multiple frames represented in a single frame block (frmb)!')
    else:
        return list(frm_cats)[0]
