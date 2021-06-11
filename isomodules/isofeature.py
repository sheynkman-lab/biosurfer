# module with class of features


class Feature():
    """Abstract class that holds methods in common across feat_objs.

       Features come in two classes. They can be either group-dependent, or
       group-independent. A group is a container of 2 or more ORFs.
    """
    def __init__(self, iso_obj, ftype):
        # self.grp -> property, grab from ref. location in featf object
        self.obj = iso_obj
        self.update_ref_to_this_feat_in_obj(ftype)

    @property
    def grp(self):
        return self.featf.grp

    def __repr__(self):
        return self.name

    def update_ref_to_this_feat_in_obj(self, ftype):
        """Add references to this newly created feat_obj to iso_obj.
           Expected input - dom, ptm, isr, etc.
           ORF attribute is orf.doms, orf.ptms, CDS attribute is cds.doms, cds.ptms, etc.
           When orf.dom requested, retrieves the dom in orf.doms matching curr.
           feat.
        """
        ftype_plural = ftype + 's'
        getattr(self.obj, ftype_plural).add(self)

class FeatureFull(Feature):
    """Represents a feature that is mapped to an ORF.
       Holds connected series of feature blocks, subblocks, residues.
       Note - sometimes blocks may not exist, depending on feat. type

       featbs - up-to-downstream feature blocks (range of same-cat, e.g. ins)
       featsbs - up-to-downstream feature subblocks (cds/exon-range-resolved)
       featrs - up-to-downstream feature residues (res-to-res)
       featps - up-to-downstream feature positions (featp-to-pos)

       Note - All feat_obj can point to a grp that holds at least two orfs.
    """
    def __init__(self, orf, ftype, desc=None):
        self.featf = self  # need for retrieving grp, for grp-dep. feat
        self.ftype = ftype  # e.g., dom, frm, isr
        self.desc = desc # more information about the feat. type
        self.blocks = []  # optional
        self.subblocks = []
        self.chain = []
        # self.first -> property
        # self.last -> property
        # self.orf -> property
        # self.full -> defined in subclass
        # self.name -> defined in subclass
        # self.grp -> defined in subclass if feat is grp-dependent
        Feature.__init__(self, orf, ftype)

    @property
    def orf(self):
        return self.obj

    @property
    def first(self):
        return self.chain[0]

    @property
    def last(self):
        return self.chain[-1]

    @property
    def name(self):
        str = 'ORF:{} type:{} {}-{}'.format(self.orf.name, self.ftype, self.first.idx, self.last.idx)
        return str


class FeatureBlock(Feature):
    """Represents a range across one or more cds/exon with same-cat feat.
       Note - this may not exist for some feature types (e.g. domain)
       Note - this represents an intermediate layer with no direct link to orf
              this info. is obtained through featf
    """
    def __init__(self, featf, featr_chain, cat):
        self.cat = cat
        self.featf = featf
        # self.first -> property
        # self.last -> property
        self.subblocks = []
        self.chain = featr_chain
        # self.ord -> property
        link_up_parent_child_references(self)

    @property
    def first(self):
        return self.chain[0]

    @property
    def last(self):
        return self.chain[-1]

    @property
    def ord(self):
        return self.featf.blocks.index(self) + 1

    @property
    def orf(self):
        return self.featf.orf

    @property
    def name(self):
        # need a name to compute ordinal on demand
        return 'frmb: {} {} -- {}'.format(self.featf.name, self.first.name, self.last.name)

    def link_up_parent_child_references(self):
        """Set upper (frmf) and lower (frmr) references."""
        # TODO: - check if correct code, b/c transferred from isocreatefeat.py
        self.frmf.blocks.append(self)
        for frmr in self.chain:
            frmr.featb = self



class FeatureSubblock(Feature):
    """A segment of a feature block, irreducable based on CDS structure.
       Note - some features lack featb, so featsb maps to featf directly

       iso_obj is cds, working with protein-centric features for now
    """
    def __init__(self, featf, cds, featr_chain):
        # self.cat -> property
        # self.ord -> property
        self.featf = featf
        self.featb = None  # update with mother feature block, if exists
        # self.first -> property
        # self.last -> property
        self.cds = cds
        # self.cds_first -> property
        # self.cds_last -> property
        self.chain = featr_chain
        # self.name -> property
        self.link_up_parent_child_references()
        Feature.__init__(self, cds, self.featf.ftype)

    @property
    def cat(self):
        if self.featb:
            return self.featb.cat
        else:
            return self.featf.cat

    @property
    def ord(self):
        return self.featf.subblocks.index(self) + 1

    @property
    def first(self):
        return self.chain[0]

    @property
    def last(self):
        return self.chain[-1]

    @property
    def name(self):
        return '{} {}-{}'.format(self.cds, self.first.res.idx, self.last.res.idx)

    def link_up_parent_child_references(self):
        # TODO: - check that function works, moved from isocreatefeat.py to here
        # on 200503
        """Update upwer/lower references."""
        self.featf.subblocks.append(self)
        for featr in self.chain:
            featr.featsb = self


class FeatureResidue(Feature):
    """A feature mapped to a residue."""
    def __init__(self, featf, res):
        # self.idx = idx  -> property, 1-base index of feature
        self.featf = featf
        self.res = res
        self.featb = None  # updated later, if exists
        self.featsb = None  # updated later, if exists
        # self.name -> defined in subclass
        Feature.__init__(self, res, self.featf.ftype)

    @property
    def idx(self):
        return self.featf.chain.index(self) + 1

    @property
    def cat(self):
        if self._cat:
            return self._cat
        else:
            return self.featf.cat

    @property
    def aa(self):
        return self.res.aa

    @property
    def name(self):
        return str(self.res.idx) + '-' + self.res.aa






class DomainFull(FeatureFull):
    """Represents a domain mapped to an ORF"""
    def __init__(self, orf, cat, acc, desc, eval, stat='direct'):
        self.cat = cat  # categories of feature (e.g. dom: dbd, reg)
        self.acc = acc  # e.g., pfam accession, linear motif acc, ptm accession
        self.desc = desc  # e.g., pfam name
        self.eval = float(eval) # -1 if non-existent
        self.stat = stat # direct (mapped domain) or transferred (from aln_obj)
        FeatureFull.__init__(self, orf, 'dom')

    @property
    def name(self):
        return self.obj.name + '|' + self.ftype + '-' + self.desc

    @property
    def full(self):
        """Char. representation of domain, cds, and AA seq. tracks."""
        chain = self.orf.res_chain
        self.orf.current_feat = self
        domain_chain = ''.join(['X' if res.dom else '-' for res in chain])
        cds_chain = ''.join(['|' if res.is_at_cds_edge else str(res.cds.ord) for res in chain])
        seq_chain = ''.join([res.aa for res in chain])
        ostr = '{:16s}{}\n{:16s}{}\n{:16s}{}'.format(self.desc, domain_chain,
                                                     'CDS ord.', cds_chain,
                                                     'AA sequence', seq_chain)
        return ostr

    def return_aln_line(self, grp, res):
        """A single 'track' showing the feat mapped to orf AA chain."""
        self.orf.current_feat = self
        chain = grp.alnf.chain
        feat_str = ''.join(['X' if getattr(aln, res).dom else '-' for aln in chain])
        self.orf.current_feat = None
        return '{:16s}{}\n'.format(self.desc, feat_str)


class DomainSubblock(FeatureSubblock):
    """A domain-specific feature subblock (domsb)."""
    def __init__(self, featf, cds, featrs):
        # domain subblock
        FeatureSubblock.__init__(self, featf, cds, featrs)

class DomainResidue(FeatureResidue):
    """A residue that represents an AA that is part of a domain."""
    def __init__(self, featf, res):
        # domain residue
        # self.orf -> property
        # self.idx -> property (1-based index of domr)
        # self.name -> property
        FeatureResidue.__init__(self, featf, res)

    @property
    def orf(self):
        return self.featf.orf

    @property
    def idx(self):
        # if issues with slowness, do-pre caching like in CDS.ord
        return self.featf.chain.index(self) + 1

    @property
    def name(self):
        return str(self.idx) + '-' + self.res.aa






class FrameFull(FeatureFull):
    """Represents the relative frame for an orf from alignment of two orfs.
       Note - assumes that this object is created immeidately after the
              creation of an alignment object, where the res.rfrm is populated.
    """
    def __init__(self, orf, grp):
        self.grp = grp
        # self.name -> property
        # self.full -> property
        FeatureFull.__init__(self, orf, 'frm')

    @property
    def name(self):
        return 'frmf: {} ({})'.format(self.obj.name, self.grp.repr_orf.name)

    @property
    def full(self):
        aa_str = ''.join([res.aa for res in self.orf.res_chain])
        self.orf.current_feat = self
        frm_str = ''.join([frmr.cat for frmr in self.chain]) # frame (1,2,3)
        stat_str = ''.join([frmr.status for frmr in self.chain]) # inferred stat.
        return aa_str + '\n' + frm_str + '\n'  + stat_str
        # TODO: - code a string representation of a frame object
        # options - frame for every AA in orf_obj
        #           coordinate range for frame block

class FrameBlock(FeatureBlock):
    """Represents a block of same-frame frmrs."""
    def __init__(self, frmf, frmr_chain, cat):
        # cat - 1, 2, or 3 (represents relative frame)
        FeatureBlock.__init__(self, frmf, frmr_chain, cat)

class FrameResidue(FeatureResidue):
    """A residue that represents the relative frame of translation.

       self.cat - direct, if the relative frame of the residue is inferred
                  indirect, if residue only in alt. orf, so inferred frame
    """
    def __init__(self, frmf, res, cat, status):
        self.cat = cat # the relative frm compared to ref. orf, values 1, 2, or 3
        self.status = self.get_frame_status(status) # whether frame is directly or indirectly inferred
        FeatureResidue.__init__(self, frmf, res)

    @property
    def name(self):
        # F1 = transposed frame 1, f1 = inferred frame 1
        frame_char = 'F'
        if self.cat == 'indirect': frame_char = 'f'
        return str(self.res.idx) + '-' + self.res.aa + ' ' + frame_char + self.cat

    def get_frame_status(self, status):
        # F = tranposed frame, f = inferred frame
        if status == '0': return 'F'
        if status == '1': return 'f'
        raise Warning('invalid frame reasidue status, needs to be 0 or 1')

# TODO: - write code for other feature types, such as ISRs

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
#         #TODO: - combine common elements of ORF.set_domain_exon_mappings and ORF_Group.set_isr_exon_mappings
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
