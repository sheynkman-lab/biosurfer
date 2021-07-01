#!/usr/bin/env python
# title           :isogroup.py
# description     :Classes that represent groupings of iso. objects.
# author          :Gloria Sheynkman
# date            :May 15th, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================


class Group():
    """General container that holds one or more iso-objects."""
    def __init__(self, objs, repr_obj=None):
        self.objs = objs
        # self.cat -> defined in subclass
        # self.name -> property
        # self.full -> full str representation
        self.update_ref_to_this_grp_in_obj()

    @property
    def name(self):
        return self.cat + ' grp: ' + '|'.join([obj.name for obj in self.objs])

    @property
    def full(self):
        # TODO: - fill in
        pass

    def update_ref_to_this_grp_in_obj(self):
        """Add reference of this newly created group to all obj.grps attr."""
        for obj in self.objs:
            obj.grps.add(self)

    def __repr__(self):
        return self.name

    def __iter__(self):
        for obj in self.objs:
            yield obj


class ORFGroup(Group):
    """Group of ORFs of arbitrary relationship."""
    def __init__(self, orfs):
        # self.orfs -> property
        # self.repr_orf -> defined in baseclass
        # self.other_orfs -> property
        # self.other_orf -> property
        Group.__init__(self, orfs)

    @property
    def orfs(self):
        return self.objs

    @property
    def other_orfs(self):
        # TODO: - test
        return [orf for orf in self.orfs if orf is not self.repr_orf]

    def enter(self):
        """Set up conditions to 'enter' or 'activate' a group."""
        for orf in self.orfs:
            orf.current_grp = self

    def exit(self):
        """Set up conditions to 'exit' or 'deactivate' a group."""
        for orf in self.orfs:
            orf.current_grp = None


class PairwiseAlignmentGroup(ORFGroup):
    """Represents a protein-sequence-centric alignment of two ORFs.

       Note - with all alignments, need to define an 'anchor' orf.
    """
    def __init__(self, anchor_orf, other_orf):
        self.cat = 'AL'
        self.repr_orf = anchor_orf
        self.anchor_orf = anchor_orf
        self.other_orf = other_orf  # consider ref to superclass ref.
        # self.anchor_orf -> property
        self.alnf = None  # update during creation 'full alignment' object
        self.frmf = None # frmf for orf2 (other_orf)
        ORFGroup.__init__(self, set([anchor_orf, other_orf]))


class SameSequenceGroup(ORFGroup):
    """Group of ORFs with the same RNA and/or protein sequence."""
    def __init__(self, orfs):
        self.cat = 'SQ'
        ORFGroup.__init__(self, orfs)


class SameIsoformStructureGroup(ORFGroup):
    """A group of ORFs with the same isoform exon-intron structure."""
    def __init__(self, orfs):
        self.cat = 'ST'
        ORFGroup.__init__(self, orfs)


class SameSequenceGencodeGroup(SameSequenceGroup):
    """Group of Gencode ENST/ENSPs that have the same protein sequence."""
    def __init__(self, orfs):
        self.cat = 'GC'
        # self.repr_orf -> property
        # self.appris_orf -> property, syn. of repr_orf
        ORFGroup.__init__(self, orfs)

    @property
    def repr_orf(self):
        """The best ORF in terms of appris and other flags."""
        flags_to_order = []  # list of [<list_of_flags>, <orf_obj>]
        for orf in self.orfs:
            flags_to_order.append([orf.flags, orf])
        repr_orf = sorted(flags_to_order)[0][1]
        return repr_orf

    @property
    def other_orfs(self):
        """Flag-ordered list of all orfs that are not the repr. orf."""
        flags_to_order = []  # list of [<list_of_flags>, <orf_obj>]
        for orf in self.orfs:
            flags_to_order.append([orf.flags, orf])
        other_orfs = [row[1] for row in sorted(flags_to_order)[1:]]
        return other_orfs

    @property
    def appris_orf(self):
        return self.repr_orf
