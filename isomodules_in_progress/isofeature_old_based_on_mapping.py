# module with class of features, which include children of all sorts of
# DNA/RNA/Protein features



class FeatureFull():
    """Represents a biophysical or functional feature mapped to an ORF.
       Called a 'full' feature or featf in attributes.
    """
    def __init__(self, orf, name):
        # self.mode -> defined in subclass (protein or rna)
        # self.cat -> defined in subclass (e.g. domain)
        self.orf = orf
        self.name = name  # e.g. domain name
        # self.first ->  property
        # self.last -> property
        self.featf = self  # needed for map.feat.featf calls
        self.blocks = []  # feat_blocks
        # self.featbs -> property, feature blocks
        self.chain = []  # chain of feat_pos or feat_res
        # self.featrs -> property, feature residues
        # self.len -> property
        self.map = None  # mapping to orf
        self.annot = None  # map to annot_obj

    @property
    def first(self):
        """First (most upstream) feat_res or feat_pos in the feat_obj."""
        return self.chain[0]

    @property
    def last(self):
        """Last (most downstream) feat_res or feat_pos in the feat_obj."""
        # TODO - return None if feat of length 1
        return self.chain[-1]

    @property
    def featbs(self):
        """List of feature blocks of this feature."""
        return self.blocks

    @property
    def featrs(self):
        """List of feature residues of this feature."""
        return self.chain

    @property
    def len(self):
        return len(self.chain)

    @property
    def full(self):
        return self.orf.name + ' -- ' + self.name

    def __repr__(self):
        return self.name

    def __iter__(self):
        for block in self.blocks:
            yield block


class Domain(FeatureFull):
    """Represents a protein domain (e.g., pfam, interproscan)."""
    def __init__(self, orf, name, acc, subcat=''):
        self.mode = 'protein'
        self.cat = 'domain'
        self.subcat = subcat  # e.g., DBD, Act, Repr
        self.acc = acc  # e.g. pfam accession
        FeatureFull.__init__(self, orf, name)


class FeatureBlock():
    """A sub-segment of a feat_obj, which maps to an exon or cds.
       Called a 'block' feature or featb.

       Note - 'pt' is either res (if protein-centric) or pos (if nuc-centric)
       Assumes that feat_pts is a consecutive chain of res or pos mapping
       to the same cds or exon.
    """
    def __init__(self, featf, feat_pts):
        self.featf = featf
        # self.cat -> property
        # self.orf -> property
        # self.ord -> property
        # self.name -> property
        self.chain = feat_pts
        # self.featrs -> property
        # TODO - map is missing!
        # self.map

    @property
    def cat(self):
        return self.featf.cat

    @property
    def orf(self):
        return self.featf.orf

    @property
    def ord(self):
        """The ordinal of the feat_block, from order of block in feat."""
        return self.featf.blocks.index(self) + 1

    @property
    def name(self):
        return self.featf.name + '|F' + str(self.ord)

    @property
    def featrs(self):
        """Return a list of feat_res_objs corr. to this block. Note - not valid
           for when the block holds feat_pos_objs.
        """
        return self.chain

    def __repr__(self):
        return self.name

    def __iter__(self):
        for pt in self.chain:
            yield pt


class FeatureResidue():
    """Represents a residue-specific point of a feat_obj."""
    def __init__(self, featf):
        # self.name -> property
        self.featf = featf
        self.featb = None  # populated upon creation of featb (featblock)
        # self.idx -> property
        self.map = map  # points to res_obj on orf

    @property
    def name(self):
        return str(self.idx) + '-' + self.map.res.aa

    @property
    def idx(self):
        return self.featf.chain.index(self) + 1

    def __repr__(self):
        return str(self.idx) + '-' + self.map.obj.aa


class FrameResidue(FeatureResidue):
    """Represents a frame of translation, at the residue-level.

       Note - determination of frame is bottom-up (frmr->frmb->frmf)
       type - 'direct' or 'inferred'
    """
    def __init__(self, frame_ord, type='direct'):
        self.type = type  # direct (pos/res-proprogated) or inferred frame
        self.ord = frame_ord
        self.map = None  # set after this feat. creation
        FeatureResidue.__init__(self, None)

    def __repr__(self):
        return self.ord


class FeaturePosition():
    """Represents a position-specific segment of a feat_obj."""
    pass


#
# class ISR(Feature):
#     """Isoform-specific region.  Defined when ogrp created."""
#     def __init__(self, orf, name, start, end):
#         Feature.__init__(self, orf, name, start, end)
#         self.cat = 'isr'
#         self.maps = [] # exon-isr mappings
#
#
# class LM(Feature):
#     # TODO- fil
#     pass
#
#
#
#
#
# # ******************
# # fractional splicing (alt. vs. const. regions)
#     def create_and_link_frac_splice_obj_to_pos(self):
#         """For each pos_obj in orf of ogrp, determine fractional splice properties,
#         create fs_obj (fractional splice) and link to pos.fs attribute."""
#         # first, make pos frequency counter, based on pos (abs) coord position
#         # so that same-coord pos can be properly grouped
#         all_pos = []
#         for orf in self.orfs:
#             for pos in orf.positions:
#                 all_pos.append(pos.coord)
#         pos_cts = Counter(all_pos)
#         num_orfs = len(self.orfs) # to calculate total number of isoforms
#         # init Frac_Splice object
#         for orf in self.orfs:
#             for pos in orf.positions:
#                 inc_ct = pos_cts[pos.coord]
#                 tot_ct = num_orfs
#                 fs_obj = Frac_Splice(inc_ct, tot_ct)
#                 pos.fs = fs_obj
#                 pos.cat = fs_obj.cat
#
#     def set_exon_frac_splice_cats(self):
#         """Assuming pos_objs are populated, set exon_obj splice category."""
#         for orf in self.orfs:
#             for exon in orf.exons:
#                 cat_chain = '' # for exon, chain of frac_splice cats (e.g. AAAAACCCCC), A=alt, S=subset, C=const.
#                 for pos in exon:
#                     cat_chain += pos.fs.cat
#                 # determine the exon category
#                 # A = alternative, M = mixed, C = constitutive
#                 if len(set(cat_chain)) == 1:
#                     repr_char = cat_chain[0]
#                     if repr_char == 'A':
#                         exon.cat = 'A'
#                     elif repr_char == 'S':
#                         exon.cat = 'M'
#                     elif repr_char == 'C':
#                         exon.cat = 'C'
#                     else:
#                         raise UserWarning('Error setting exon frac_splice category:' + exon.name)
#                 else:
#                     exon.cat = 'M'
#
# class Frac_Splice():
#     """Object holds information related to the isoform-specificity of a pos_obj."""
#     def __init__(self, inc, all):
#         self.inc = inc # number of orfs in ogrp that pos is represented in
#         self.all = all # total number of orfs a pos can be in
#         self.frac = inc/float(all)
#         self.frac_abund = None # abundance-corrected version of fs.frac, populate later
#         self.cat = self.set_fs_cat()
#
#     def __repr__(self):
#         ostr = '{:^2} / {:^2} = {:^3.1f}, {:^6}'.format(self.inc, self.all, self.frac, self.cat)
#         return ostr
#
#     def set_fs_cat(self):
#         """Set fs category based on an enum.
#         Category can be:
#         Cat.A (alternative, pos is specific to one ORF)
#         Cat.S (subset, pos is in subset of ORFs)
#         Cat.C (constitutive, pos is in all ORFs of ogrp)
#         """
#         if self.frac == 1:
#             return 'C'
#         elif self.inc == 1:
#             return 'A'
#         else:
#             return 'S'
#
# # ********************
# # isoform-specific regions (ISRs)
#
#     def define_create_and_map_isrs(self):
#         """Find ISRs of ORFs in ogrp.
#         ISR = isoform-specific region, defined as run of pos.cat=A (alternative)
#         ISR based on orf.idx range, then create ISR object, then ISR-exons mapping objects.
#         Link ISR-exons via mapping object."""
#
#         # methods used in this func.
#         def get_cat_chain(orf):
#             """Return a string repr. of the alt/sub/con categories for ordered pos_obj in orf.
#             e.g. 'CCCCCCCAAAAAAASSSSSCCCCAA'"""
#             cat_chain = ''
#             for pos in orf:
#                 cat_chain += pos.cat
#             return cat_chain
#
#         def explode_cat_chain(cat_chain):
#             """Split cat_chain base on runs of same-char."""
#             return [item[0] for item in re.findall(r"((.)\2*)", cat_chain)]
#
#         def get_isr_ranges(orf, blocks):
#             """From list of blocks, return dict of orf 1-based idx ranges of all alt blocks.
#             e.g. ['CCC', 'AA', 'S', 'AAA'] returns {'ORF1-ISR01':[4,5], 'ORF1-ISR02':[7,9]}"""
#             isr_ranges = {}
#             st_idx = 1 # start 1-based orf idx
#             en_idx = None # end 1-baed orf idx
#             isr_idx = 1
#             for block in blocks:
#                 blen = len(block)
#                 en_idx = st_idx + blen - 1
#                 if 'A' in block:
#                     acc = '{}-ISR{:02d}'.format(orf.name, isr_idx)
#                     isr_idx += 1
#                     isr_ranges[acc] = [st_idx, en_idx]
#                 st_idx = en_idx + 1
#             return isr_ranges
#
#         #TODO - combine common elements of ORF.set_domain_exon_mappings and ORF_Group.set_isr_exon_mappings
#         def set_isr_exon_mappings(orf):
#             """For each isr, encode isr-exon mappings through a Mapping object.
#             Note: code similar to set_domain_exon_mappings under ORF class."""
#             for isr in orf.isrs:
#                 isr_cur_idx = 1 # track isr ranges (as split by exons)
#                 for exon in orf.exons:
#                     exon_start = 0
#                     exon_end = 0
#                     for pos in exon:
#                         if isr.start < pos.idx < isr.end:
#                             if pos.ridx == 1:
#                                 exon_start = pos.ridx
#                             if pos.ridx == exon.len:
#                                 exon_end = exon.len
#                         if pos.idx == isr.start:
#                             exon_start = pos.ridx
#                         if pos.idx <= isr.end:
#                             exon_end = pos.ridx
#                     if exon_start and exon_end:
#                         sub_isr_len = exon_end - exon_start
#                         isr_start = isr_cur_idx
#                         isr_end = isr_cur_idx + sub_isr_len
#                         isr_cur_idx = isr_end + 1
#                         isr_dom_mapping = Mapping(exon, isr, exon_start, exon_end, isr_start, isr_end)
#                         isr.maps.append(isr_dom_mapping)
#                         exon.maps.append(isr_dom_mapping)
#                     # else:
#                     #     raise UserWarning('error setting isr mappings:' + exon.name + ',' + isr.name)
#
#         for orf in self.orfs:
#             chain = get_cat_chain(orf) # e.g. 'CCCAASS'
#             blocks = explode_cat_chain(chain) # e.g. ['CCC', 'AA', 'SS']
#             isr_dict = get_isr_ranges(orf, blocks)
#             for acc, pair in isr_dict.items():
#                 st, en = pair # start/end of isr range (1-base of orf)
#                 isr_obj = ISR(orf, acc, st, en)
#                 orf.isrs.append(isr_obj)
#                 set_isr_exon_mappings(orf)
