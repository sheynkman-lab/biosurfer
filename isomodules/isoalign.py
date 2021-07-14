#!/usr/bin/env python
# title           :isoalign.py
# description     :Classes that represent alignments between iso. objects.
#                  Also, functions to compute and print-out aa splice aligns.
# author          :Gloria Sheynkman
# date            :May 22nd, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================


class Alignment():
    """Abstract class that holds methods in common across aln_objs."""
    def __init__(self, anchor_obj, other_obj):
        # self.grp -> property, grab from ref. location in alnf object
        self.anchor_obj = anchor_obj
        self.other_obj = other_obj
        self.update_ref_to_this_aln_in_objs()

    # @property
    # def grp(self):
    #     return self.alnf.grp
    
    # @grp.setter
    # def grp(self, grp):
    #     self.alnf.grp = grp

    def __repr__(self):
        # previously, had the isoalign obj. name
        # class_name = self.__class__.__name__
        anchor = self.anchor_obj.name
        other = self.other_obj.name
        return '{}|{}'.format(anchor, other)

    def update_ref_to_this_aln_in_objs(self):
        """Add references to this newly created aln_obj to both objs.
           Note - not done for AlignmentBlock, which has no pt. to obj.
        """
        if self.anchor_obj:
            self.anchor_obj.alns.add(self)
        if self.other_obj:
            self.other_obj.alns.add(self)


class AlignmentFull(Alignment):
    """Represents a series of alignments across an ORF.
       Holds connected series of alignment blocks, subblocks, residues.

       alnbs - up-to-downstream alignment blocks (range of same-cat, e.g. ins)
       alnsbs - up-to-downstream alignment subblocks (cds-range-resolved)
       alnrs - up-to-downstream alignment residues (res-to-res)

       Note - In all alignments, there is a concept of the 'anchor' orf and
              the 'other' orf. Identity of anchor/other dictates aln type.

       Note - All aln_obj must point to a grp that holds two orfs.
    """
    def __init__(self, grp, anchor_orf, other_orf):
        self.alnf = self  # need for retrieving grp
        self.grp = grp
        # self.anchor_orf -> property
        # self.orf1 -> property, syn. of anchor_orf
        # self.other_orf -> property
        # self.orf2 -> property, syn. of other_orf
        self.protblocks = []  # protein-change blocks
        self.blocks = []  # nt-change blocks
        self.subblocks = []
        self.chain = []  # alnr_obj added during creation
        # self.seq1 -> property, sequence of AAs
        # self.seq2 -> property, sequence of AAs
        # self.full -> property
        # self.match_block_chain -> property
        Alignment.__init__(self, anchor_orf, other_orf)

    @property
    def anchor_orf(self):
        return self.anchor_obj

    @property
    def orf1(self):
        return self.anchor_orf

    @property
    def other_orf(self):
        return self.other_obj

    @property
    def orf2(self):
        return self.other_orf

    @property
    def seq1(self):
        return ''.join([alnr.res1.aa for alnr in self.chain])

    @property
    def seq2(self):
        return ''.join([alnr.res2.aa for alnr in self.chain])

    @property
    def cds1(self):
        return ''.join(['|' if alnr.res1.is_at_cds_edge else str(alnr.res1.cds.ord)[0] for alnr in self.chain])

    @property
    def cds2(self):
        return ''.join(['|' if alnr.res2.is_at_cds_edge else str(alnr.res2.cds.ord)[0] for alnr in self.chain])


    @property
    def full(self):
        """String representation of the orf-orf splice-based alignment."""
        cds1 = ''.join([str(alnr.res1.cds.ord)[-1] for alnr in self.chain])
        aa1 = ''.join([alnr.res1.aa for alnr in self.chain])
        aa2 = ''.join([alnr.res2.aa for alnr in self.chain])
        cds2 = ''.join([str(alnr.res2.cds.ord)[-1] for alnr in self.chain])
        frm2 = ''.join([alnr.res2.rfrm for alnr in self.chain])
        alnp_chain = ''.join([alnr.alnpb.cat for alnr in self.chain])
        alnb_chain = ''.join([alnr.alnb.cat for alnr in self.chain])
        ostr = ('{gene} {strand}\n{orf1:16s} AA:{aa1}\n'
                '               CDS1:{cds1}\n'
                '{orf2:16s} AA:{aa2}\n'
                '               CDS2:{cds2}\n'
                '                 FM:{frm2}\n'
                '               BCAT:{aln}\n'
                '               PCAT:{alnp}')
        ostr = ostr.format(gene=self.orf1.gene.name, strand=self.orf1.strand,
                        orf1=self.orf1.name, aa1=aa1, orf2=self.orf2.name,
                        aa2=aa2, frm2=frm2, aln=alnb_chain, alnp=alnp_chain,
                        cds1=cds1, cds2=cds2)
        return ostr

    @property
    def match_block_chain(self):
        """String repr. of the match-types.
           e.g. MMMMMDDDDDMMMM, returns MDM
        """
        return ''.join([alnb.cat for alnb in self.blocks])

    def get_subset_aln_chain(self, first_alnr, last_alnr):
        """Return a list of alnr between first/last_alnr.
           Note - assumes that both alnr are in self.chain and first < last.
        """
        subchain = []
        in_subset = False
        for alnr in self.chain:
            if alnr == first_alnr:
                subchain.append(alnr)
                in_subset = True
            elif alnr == last_alnr:
                subchain.append(alnr)
                in_subset = False
            elif in_subset:
                subchain.append(alnr)
            else:
                pass
        return subchain



class AlignmentProteinBlock(Alignment):
    """Represents a protein-centric change between isoforms.

       Cat - M=match, I=insertion, D=deletion, S=substitution
    """
    def __init__(self, cat, alnf, alnbs):
        self.cat = cat  # aln type, see notes above
        self.alnf = alnf  # ref. to the full align (length of orf)
        self.grp = self.alnf.grp
        # self.first -> property
        # self.last -> property
        self.blocks = alnbs
        # self.subblocks -> property
        # self.chain -> property
        # self.res_chain1 -> property
        # self.aa1 -> property
        # self.res_chain2 -> property
        # self.aa2 -> property
        self.update_references_to_parent_and_child_objects()
        Alignment.__init__(self, None, None)  # no direct anchor/other object

    @property
    def first(self):
        return self.chain[0]

    @property
    def last(self):
        return self.chain[-1]

    @property
    def subblocks(self):
        # TODO: - check
        alnsbs = []
        for alnb in self.blocks:
            for alnsb in alnb.subblocks:
                if alnsb not in alnsbs:
                    alnsbs.append(alnsb)
        return alnsbs

    @property
    def chain(self):
        alnrs = []
        for alnb in self.blocks:
            for alnr in alnb.chain:
                alnrs.append(alnr)
        return alnrs

    @property
    def res_chain1(self):
        return [alnr.res1 for alnr in self.chain if alnr.res1.aa != '-']

    @property
    def aa1(self):
        return ''.join([res.aa for res in self.res_chain1])

    @property
    def res_chain2(self):
        return [alnr.res2 for alnr in self.chain if alnr.res2.aa != '-']

    @property
    def aa2(self):
        return ''.join([res.aa for res in self.res_chain2])

    def __repr__(self):
        orf1_name = '-' if not self.first.res1 else self.first.res1.orf.name
        orf2_name = '-' if not self.first.res2 else self.first.res2.orf.name
        return 'alnpb: {}|{} {} {}-{}'.format(orf1_name, orf2_name, self.cat,
                                                    self.first.idx,
                                                    self.last.idx)

    def update_references_to_parent_and_child_objects(self):
        # set lower rand upper references
        # TODO: - check if works, transferred from isocreatealign
        self.alnf.protblocks.append(self)
        for alnb in self.blocks:
            alnb.alnpb = self
        for alnsb in self.subblocks:
            alnsb.alnpb = self
        for alnr in self.chain:
            alnr.alnpb = self




class AlignmentBlock(Alignment):
    """Represents a range across one or more cds/exon with same-cat align.

       Cat - M=match, I=insertion, D=deletion, F=frameshift
    """
    def __init__(self, cat, alnf, alnr_chain):
        self.cat = cat  # aln type, see notes above
        self.alnf = alnf  # ref. to the full align (length of orf)
        self.alnpb = None  # assigned during creation
        # self.first -> property
        # self.last -> property
        self.subblocks = []
        self.chain = alnr_chain
        self.update_references_to_parent_and_child_objects()
        Alignment.__init__(self, None, None)  # no direct anchor/other object

    @property
    def first(self):
        return self.chain[0]

    @property
    def last(self):
        return self.chain[-1]

    def __repr__(self):
        orf1_name = '-' if not self.first.res1 else self.first.res1.orf.name
        orf2_name = '-' if not self.first.res2 else self.first.res2.orf.name
        return 'alnb: {}|{} {} {}-{}'.format(orf1_name, orf2_name, self.cat,
                                                    self.first.idx,
                                                    self.last.idx)

    @property
    def ord(self):
        return self.alnf.blocks.index(self) + 1

    def update_references_to_parent_and_child_objects(self):
        # set lower and upper references
        # TODO: - check if works, transferred from isocreatealign
        self.alnf.blocks.append(self)
        for alnr in self.chain:
            alnr.alnb = self



class AlignmentSubblock(Alignment):
    """A segment of an alignment block, irreducable based on CDS structure."""
    def __init__(self, alnf, anchor_cds, other_cds, alnr_chain):
        # self.cat -> property
        # self.ord -> property
        # self.anchor_cds -> property
        # self.cds1 -> property, syn. of anchor_cds
        # self.other_cds -> property
        # self.cds2 -> property, syn. of other_cds
        self.alnf = alnf
        self.grp = self.alnf.grp
        self.alnpb = None  # updated upon isoalign creation
        self.alnb = None  # update with mother alignment block
        # self.first -> property
        # self.last -> property
        self.chain = alnr_chain
        # self.alnrs -> property, syn of chain
        self.update_references_to_parent_and_child_objects()
        Alignment.__init__(self, anchor_cds, other_cds)

    @property
    def cat(self):
        return self.alnb.cat

    @property
    def ord(self):
        return self.alnf.subblocks.index(self) + 1

    @property
    def anchor_cds(self):
        return self.anchor_obj

    @property
    def cds1(self):
        return self.anchor_cds

    @property
    def other_cds(self):
        return self.other_obj

    @property
    def cds2(self):
        return self.other_cds

    @property
    def first(self):
        return self.chain[0]

    @property
    def last(self):
        return self.chain[-1]

    @property
    def alnrs(self):
        return self.chain

    def __repr__(self):
        orf1_name = '-' if not self.first.res1 else self.first.res1.orf.name
        orf2_name = '-' if not self.first.res2 else self.first.res2.orf.name
        return 'alnsb: {}|{} {} {}-{}'.format(orf1_name, orf2_name, self.cat,
                                              self.first.idx, self.last.idx)

    def update_references_to_parent_and_child_objects(self):
        """Grab the alnb associated with this alnsb and link up refs."""
        # TODO: - check for correctness
        alnb = self.get_alnb_mapped_to_alnr_chain(self.chain)
        self.alnb = alnb
        alnb.subblocks.append(self)
        self.alnf.subblocks.append(self)
        for alnr in self.chain:
            alnr.alnsb = self


    def get_alnb_mapped_to_alnr_chain(self, alnr_chain):
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




class AlignmentResidue(Alignment):
    """An alignment between two residues."""
    def __init__(self, alnf, cat, anchor_res, other_res):
        self.cat = cat  # residue-level alignment data
        # self.idx = idx  -> property, 1-base index of alignment (NOT residue)
        # self.anchor_res -> property
        # self.res1 -> property, syn of anchor_res
        # self.other_res -> property
        # self.res2 -> property, syn of other_res
        self.alnf = alnf
        self.grp = self.alnf.grp
        self.alnpb = None  # updated during isoalign creation
        self.alnb = None  # see above
        self.alnsb = None  # see above
        Alignment.__init__(self, anchor_res, other_res)

    @property
    def idx(self):
        # TODO: - test
        return self.alnf.chain.index(self) + 1

    @property
    def anchor_res(self):
        return self.anchor_obj

    @property
    def res1(self):
        return self.anchor_res

    @property
    def other_res(self):
        return self.other_obj

    @property
    def res2(self):
        return self.other_res
