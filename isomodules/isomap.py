# This module likely will be deprecated, because I decided to directly
# map the feature object types to the ORF objects, rather than going through
# a mapping object. This is because adding the mapping object seemed to
# complicated the code.
# 200503

# module with class of mapping type, which links molec. obj with features
# or annotations


class Mapping():
    """Represents a mapping between a biomolecular object (i.e., iso_obj) and a
       feature object (feat_obj). Superclass for different types of mappings.
       Assumes that obj-feat alignment is 'monotonic', or linear and concord.

       Upon instantiation, the 'maps' attr in the respective iso_obj
       and feat_obj are updated with the newly created Mapping object, thereby
       completing the obj-map-feat linkage.

       Note - cat not defined because it's implicit which two obj. are linked
    """
    def __init__(self, obj, feat, grp=None):
        self.obj = obj  # has synonyms in subclasses
        self.feat = feat
        # self.type -> defined in subclass (range or point)
        # self.mode -> defined in subclass (protein or nucleotide)
        self.grp = grp  # set to the grp_obj for grp-specific mappings
        self.update_ref_to_this_map_in_obj_and_feat()

    def __repr__(self):
        return 'Map:' + self.obj.name + ' -- ' + self.feat.name

    def update_ref_to_this_map_in_obj_and_feat(self):
        """Add references of this newly created map_obj to obj and feat.
           Obj always holds a set of maps, feat always points to one map.
        """
        self.obj.maps.add(self)
        self.feat.map = self


class RangeMapping(Mapping):
    """A mapping between an iso_obj and a feat_obj, with matched ranges.
       Can be a full length mapping (e.g., orf_obj to feat_obj) or a subsegment
       mapping (e.g., cds_obj to feat_block).

       'pt' means res_obj, pos_obj, feat_res_obj, feat_pos_obj
           If the feat_obj is protein-centric, res and feat_res is used
           If the feat_obj is rna-centric, pos and feat_pos is used.
    """
    def __init__(self, obj, feat, pt_objs, feat_pt_objs, grp=None):
        self.type = 'range'
        self.obj_first = pt_objs[0]
        self.obj_last = pt_objs[-1]
        self.feat_first = feat_pt_objs[0]
        self.feat_last = feat_pt_objs[-1]
        # self.full -> defined in subclass, char. visualize mapping
        Mapping.__init__(self, obj, feat, grp)


class ORF_Featf_Mapping(RangeMapping):
    """A mapping between an orf and featf."""
    def __init__(self, orf, featf, pt_objs, feat_pt_objs, grp=None):
        self.orf = orf
        self.featf = featf
        # self.orf_start -> property
        # self.orf_end -> property
        # self.featf_start -> property
        # self.featf_end -> property
        # self.full -> property
        RangeMapping.__init__(self, orf, featf, pt_objs, feat_pt_objs, grp)

    @property
    def orf_start(self):
        return self.obj_start

    @property
    def orf_end(self):
        return self.obj_end

    @property
    def featf_start(self):
        return self.feat_start

    @property
    def featf_end(self):
        return self.feat_end

    @property
    def full(self):
        """Shows feature-track and iso_obj-track as lines of chars."""
        feat_str = '{:12s}'.format(self.featf.name)
        res_str = '{:12s}'.format(self.orf.name)
        for res in self.orf.res_chain:
            if self.featf in [map.featr.featf for map in res.maps]:
                feat_str += 'X'
            else:
                feat_str += '-'
            res_str += res.aa
        return feat_str + '\n' + res_str


class CDS_Featb_Mapping(RangeMapping):
    """A mapping between a cds and sub-segment of feat."""
    def __init__(self, cds, featb, pt_objs, feat_pt_objs, grp=None):
        self.cds = cds
        self.featb = featb  # feature 'block', i.e., cds-level feat. segment
        # self.cds_start -> property
        # self.cds_end -> property
        # self.featb_start -> property
        # self.featb_end -> property
        RangeMapping.__init__(self, cds, featb, pt_objs, feat_pt_objs, grp)

    @property
    def cds_start(self):
        return self.obj_start

    @property
    def cds_end(self):
        return self.obj_end

    @property
    def featb_start(self):
        return self.feat_start

    @property
    def featb_end(self):
        return self.feat_end

    @property
    def full(self):
        """Shows feature-track and iso_obj-track as lines of chars.
           Specific to cds_obj-featblock_obj mapping.
        """
        feat_str = '{:16s}'.format(self.featb)
        res_str = '{:16s}'.format(self.cds)
        for res in self.cds.chain_trimmed:
            if self.featb in [map.featr.featb for map in res.maps]:
                feat_str += 'X'
            else:
                feat_str += '-'
            res_str += res.aa
        return feat_str + '\n' + res_str


class PointMapping(Mapping):
    """A mapping between a pos/res object and a pos/res feat. object."""
    def __init__(self, obj_pt, feat_pt, grp=None):
        self.type = 'point'
        Mapping.__init__(self, obj_pt, feat_pt, grp)


class Res_Featr_Mapping(PointMapping):
    """A mapping between a res object and a feat_res."""
    def __init__(self, obj_pt, feat_pt, grp=None):
        self.mode = 'protein'
        self.res = obj_pt
        self.featr = feat_pt
        PointMapping.__init__(self, obj_pt, feat_pt, grp)


class Pos_Featp_Mapping(PointMapping):
    """A mapping between a pos object and a feat_pos."""
    def __init__(self, obj_pt, feat_pt, grp=None):
        self.mode = 'nucleotide'
        self.pos = obj_pt
        self.featp = feat_pt
        PointMapping.__init__(self, obj_pt, feat_pt, grp)
