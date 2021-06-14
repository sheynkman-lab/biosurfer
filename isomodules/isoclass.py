#!/usr/bin/env python
# title           :isoclass.py
# description     :Classes representing isoform-related objects.
# author          :Gloria Sheynkman
# date            :May 1st, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

from Bio.Seq import Seq
import math
from collections import Counter
from collections import defaultdict
import itertools

"""Note - three core assumptions of all isoforms per gene:
          1) are protein-codng
          2) have good alignments to the genome
          3) are of the same strand
"""

class Biomolecule():
    """Abstract superclass for all biomolecular objects.
       Place to keep shared functions and attributes.
    """
    def __init__(self):
        # self.chrom -> property
        # self.strand -> property
        # self.len -> property

        # grp-independent features, active-feature-dependent retrieval
        self.doms = set()
        self.ptms = set()
        # grp-dependent features, active-grp-dependent retrieval
        self.alns = set()
        self.frms = set()
        self.isrs = set()

    #TODO: - determine what the setter is for - had errored out when
    # calling ORF.chrom and ORF.strand when it was from Gencode Gene
    @property
    def chrom(self):
        if hasattr(self, 'gene'):
            return self.gene.chrom
        else:
            return self._chrom

    @chrom.setter
    def chrom(self, value):
        self._chrom = value

    @property
    def strand(self):
        if hasattr(self, 'gene'):
            return self.gene.strand
        else:
            return self._strand

    @strand.setter
    def strand(self, value):
        self._strand = value

    @property
    def len(self):
        return len(self.seq)

    def __repr__(self):
        return self.name

    def __len__(self):
        return len(self.seq)

    @property
    def aln(self):
        """Return a single aln_obj corresponding to the 'active' group."""
        current_grp = self.orf.current_grp
        for aln in self.alns:
            if aln.grp == current_grp:
                return aln
        return None

    @property
    def dom(self):
        """Return a single dom_obj corresponding to the 'active' dom_obj."""
        current_feat = self.orf.current_feat
        for dom in self.doms:
            if dom.featf == current_feat:
                return dom
        return None

    @property
    def ptm(self):
        """Return a single ptm_obj corresponding to the 'active' ptm_obj."""
        current_feat = self.orf.current_feat
        for ptm in self.ptms:
            if ptm.featf == current_feat:
                return ptm 
        return None

    # @property
    # def frm(self):
    #     """Return a single frm_obj corresponding to the 'active' frm_obj."""
    #     current_feat = self.orf.current_feat
    #     for frm in self.frms:
    #         if frm.featf == current_feat:
    #             return frm
    #     return None
    @property
    def frm(self):
        return self._frm

    @frm.setter
    def frm(self, frm_value):
        self._frm = frm_value





class Gene(Biomolecule):
    """Abstract class for all base classes of type Gene.
       Assume biotype is 'protein_coding'. Assume all sub-obj. same strand.
       cat, start, and end are attr required to be set by Sub class.
    """
    def __init__(self, name, chrom, strand):
        # self.cat -> defined in sub class, category of "Gene"
        self.name = name  # symbol
        self.chrom = chrom  # chr1, chr2, etc.
        self.strand = strand
        # self.start -> defined in sub class, left-most pos
        # self.end -> defined in sub class, right-most pos
        # self.first -> property, up-most pos among orfs
        # self.last -> property, dn-most pos among orfs
        self.orfs = set()
        # self.exons = set() -> property
        # self.poss = set() -> property
        # self.cdss = set() -> property
        # self.ress = set() -> property
        Biomolecule.__init__(self)

    @property
    def first(self):
        """Upstream-most pos_obj, among all orfs in gene.
           When multi. pos_obj same upstream coord, report one random.
        """
        poss = []  # list of [coord, pos_obj]
        for orf in self.orfs:
            poss.append([orf.first.coord, orf.first])
        if self.strand == '+':
            first = sorted(poss)[0][1]
        elif self.strand == '-':
            first = sorted(poss)[-1][1]
        return first

    @property
    def last(self):
        """Downstream-most pos_obj, among all orfs in gene.
           When multi. pos_obj same downstream coord, report one random.
        """
        poss = []  # list of [coord, pos_obj]
        for orf in self.orfs:
            poss.append([orf.last.coord, orf.last])
        if self.strand == '+':
            last = sorted(poss)[-1][1]
        elif self.strand == '-':
            last = sorted(poss)[0][1]
        return last

    @property
    def exons(self):
        """Return all exons assoc. with gene."""
        exons = set()
        for orf in self:
            exons.update(orf.exons)
        return exons

    @property
    def poss(self):
        """Return all pos_obj assoc. with gene."""
        poss = set()
        for orf in self:
            poss.update(orf.chain)
        return poss

    @property
    def cdss(self):
        """Return all cdss assoc. with gene."""
        cdss = set()
        for orf in self:
            cdss.update(orf.cdss)
        return cdss

    @property
    def ress(self):
        """Return all res-obj assoc. with gene."""
        ress = set()
        for orf in self:
            ress.update(orf.ress)
        return ress

    @property
    def orf_pairs(self):
        """Return all possible pairs of orfs, as a list of lists.
           e.g., [1, 2, 3] -> [1, 2], [1, 3], [2, 3]
        """
        orf_pairs = []
        for (orf1, orf2) in itertools.combinations(self.orfs, 2):
            orf_pairs.append([orf1, orf2])
        return orf_pairs

    @property
    def orf_pairs_permute(self):
        """Return all possible orf pairs, permutations, as a list of lists.
           e.g. [1, 2, 3] -> [1, 2], [2, 1], [1, 3], etc.
        """
        orf_pairs = []
        for (orf1, orf2) in itertools.permutations(self.orfs, 2):
            orf_pairs.append([orf1, orf2])
        return orf_pairs

    @property
    def same_seq_orfs(self):
        """Return groups of same-protein-sequence orf_objs.
           e.g., [[orf1_seq1, orf2_seq1], [orf3_seq2]]
        """
        seqs = defaultdict(list)
        for orf in self:
            seqs[orf.tr].append(orf)
        return seqs.values()

    @property
    def orfs_len_ordered_desc(self):
        """Return a list of orfs, order from longest to shortest."""
        orfs = sorted([[len(orf.seq), orf] for orf in self.orfs], reverse=True)
        return [x[1] for x in orfs]

    def __getitem__(self, isoname):
        for orf in self.orfs:
            if orf.name == isoname:
                return orf

    def __iter__(self):
        for orf in self.orfs:
            yield orf

    def __len__(self):
        return abs(self.start - self.end) + 1

    def full(self):
        ostr = '{} {} {}:{}-{}'.format(self.name, self.strand, self.chrom,
                                       self.start, self.end)
        return ostr

    def get_exon(self, exon_name):
        for exon in self.exons:
            if exon.name == exon_name:
                return exon

    def get_other_orfs(self, anchor_orf):
        """Return all other orfs of a gene that is not an anchor_orf."""
        other_orfs = [orf for orf in self.orfs if orf != anchor_orf]
        return other_orfs



class GencodeGene(Gene):
    """Represents a gen_obj based on a Gencode annotation.

       Note - Since start and end absolute coords. are known, assigned now.
    """
    def __init__(self, ensg, name, chrom, strand, start, end):
        self.ensg = ensg
        self.cat = 'GC'
        self.start = start  # absolute coord. on genome
        self.end = end  # see above
        # self.appris_orf -> property
        # self.repr_orf -> syn. of appris_orf
        # self.redundant_seq_orfs  # on-the-fly attr. to stash red-prot-seq orf
        Gene.__init__(self, name, chrom, strand)

    @property
    def appris_orf(self):
        """Define orf which is appris principle, based on appris tags.
           In case of ties or absense of appris annot., choose orf with highest
           ranked transcript annots (transcript_support_level, etc.).
        """
        flags_to_order = []  # list of [<list_of_flags>, <orf_obj>]
        for orf in self:
            flags_to_order.append([orf.flags, orf])
        appris_orf = sorted(flags_to_order)[0][1]
        return appris_orf

    @property
    def repr_orf(self):
        return self.appris_orf

    @property
    def other_orfs(self):
        other_orfs = set()
        for orf in self.orfs:
            if orf != self.repr_orf:
                other_orfs.add(orf)
        return other_orfs

    @property
    def ref_alt_pairs(self):
        """Return a list of gencode ref/alt orf_obj pairs. ref is best
           appris. Pairs are returned by ranking by isoname.
        """
        pairs = []  # [[appris_orf_obj, orf_obj], [appris_orf_obj, orf_obj2]..]
        for orf in self.orfs:
            if orf != self.appris_orf:
                pairs.append([self.appris_orf, orf])
        return sorted(pairs)



class SequenceBasedGene(Gene):
    """Represents a gen_obj based on independently generated full-length
       transcript sequence (e.g. full-length sequencing, isoform cloning).
    """
    def __init__(self, cat, name, chrom, strand):
        self.cat = cat  # e.g. CL, PB
        # self.start -> property
        # self.end -> property
        Gene.__init__(self, name, chrom, strand)

    @property
    def start(self):
        pass

    @property
    def end(self):
        pass


class ORF(Biomolecule):
    """Abstract class for all base classes of type ORF.

       Note - All child objects (exon, junc, cds, pos, res) in lists in upstream
       to downstream order.
    """
    def __init__(self, name, gene):
        # self.cat -> set in base class (e.g. GC, PB, CL)
        self.orf = self  # to grab current_feat/current_grp in superclass
        self.name = name
        self.gene = gene
        # self.chrom -> defined in Biomolecule
        # self.strand -> defined in Biomolecule
        # self.first -> property, upstream-most pos_obj
        # self.last -> property, downstream-most pos_obj
        self.seq = ''  # ATG-stop (can have early stop), set in pop. seq-rel. obj.
        # self.tr -> property
        # self.frac_tr -> property
        self.exons = []  # exon_obj get added during isocreate
        # self.cdss -> property
        # self.juncs -> property
        # self.chain -> property
        # self.res_chain -> property
        self.res_chain_populated = []  # if res_chain populated, call this
        self.current_grp = None  # state-specific attr, current grp orf is in
        self.current_feat = None  # state-specific attr, current feat
        self.grps = set()  # holds grp_obj where orf belongs
        self.plot_ord = 1  # optional ordinal to define plot ordering (can set it later)
        Biomolecule.__init__(self)

    @property
    def first(self):
        """Upstream-most pos_obj in orf."""
        return self.exons[0][0]

    @property
    def last(self):
        """Downstream-most pos_obj in orf."""
        return self.exons[-1][-1]

    @property
    def tr(self):
        """Translation of self.seq, ATG-to-Stop."""
        return str(Seq(self.seq).translate(to_stop=True))

    @property
    def frac_tr(self):
        """Fraction of input ORF that is translated."""
        num_translated_nt = len(self.tr) * 3
        num_total_nt = len(self.seq)
        frac_trans = float(num_translated_nt)/float(num_total_nt)
        return frac_trans

    @property
    def cdss(self):
        """List of cds_objs, upstream-to-downstream ordering."""
        cdss = [exon.cds for exon in self.exons if exon.cds]
        return cdss

    @property
    def ress(self):
        """Return all res_obj assoc. with this orf."""
        ress = set()
        for cds in self.cdss:
            ress.update(cds.chain)
        return ress

    @property
    def chain(self):
        """Gather pos_obj among exons, and return as list."""
        pos_chain = []
        for exon in self.exons:
            pos_chain.extend(exon.chain)
        return pos_chain

    @property
    def res_chain(self):
        """Gather res_obj among cdss, and return as list."""
        if self.res_chain_populated:
            return self.res_chain_populated
        else:
            res_chain = []
            for cds in self.cdss:
                res_chain.extend(cds.chain_trimmed)
            self.res_chain_populated = res_chain
            return res_chain

    @property
    def juncs(self):
        juncs = []  # up-to-dn junc_obj
        for exon in self.exons:
            if exon.dn_junc:
                juncs.append(exon.dn_junc)
        return juncs

    @property
    def cds_chain(self):
        """List of res_obj of the orf, ordered upstream-to-downstream."""
        res_chain = []
        for cds in self.cdss:
            res_chain.extend(cds.chain)
        return res_chain

    @property
    def blen(self):
        """Cumulative 'block' length of exon ranges."""
        blen = 0
        for exon in self.exons:
            blen += exon.blen
        return blen

    @property
    def full(self):
        return self.name + ' ' + self.seq

    def __iter__(self):
        for pos in self.chain:
            yield pos

    def __getitem__(self, i):
        """Fetching index 1 means getting the first position."""
        return self.chain[i-1]

    def __getslice__(self, i, j):
        """Slicing of pos out of orf is 0-based."""
        return self.__getitem__(slice(i, j))

    def __hash__(self):
        """Need to redefine hash func. if redefining __eq__."""
        return id(self)

    def __lt__(self, other):
        return self.name < other.name

    def __le__(self, other):
        return self.name <= other.name

    def __eq__(self, other):
        return self.name == other.name

    def __ne__(self, other):
        return self.name != other.name

    def __gt__(self, other):
        return self.name > other.name

    def __ge__(self, other):
        return self.name >= other.name

    def get_exon(self, exon_name):
        for exon in self.exons:
            if exon.name == exon_name:
                return exon

    def set_rel_start_and_end(self):
        """This is probably for the iso-imager."""
        self.rel_start = self.exons[0].rel_start
        self.rel_end = self.exons[-1].rel_end


class GencodeORF(ORF):
    """Represents a Genocode-based orf_obj.

       flags - annotations extracted from Gencode gtf
               used to rank orfs of a gene to find appris principal orf
               of the form -> [appris, level, tsl, rank, enst_idx]
                  appris - appris annotation (e.g. appris_principal_1)
                  level - transcript evidence level
                  tsl - transcript support level
                  rank - index from isoform name (e.g. for ZF-201, is 201)
                  enst_idx - ENST number
       is_basic - is part of Gencode basic set (versus comprehensive set)
    """
    def __init__(self, start, end, enst, ensp, appris, is_basic, flags,
                 isoname, gene, start_status, end_status):
        self.cat = 'GC'
        self.start = start
        self.end = end
        self.enst = enst
        self.ensp = ensp
        self.is_basic = is_basic
        self.cds_start_nf = start_status
        self.cds_end_nf = end_status
        self.flags = flags  # [appris, level, tsl, rank, enst_idx]
        # self.appris -> property
        # self.level -> property
        # self.tsl -> property
        ORF.__init__(self, isoname, gene)

    @property
    def appris(self):
        return self.flags[0]

    @property
    def level(self):
        return self.flags[1]

    @property
    def tsl(self):
        return self.flags[2]


class SequencedORF(ORF):
    """An experimentally-derived ORF, either from full-length sequencing or
       isoform cloning pipelines.
    """
    def __init__(self, cat, start, end, isoname, gene, orfid=None):
        self.cat = cat
        self.start = start
        self.end = end
        self.orfid = orfid  # optional, for 6K collection
        ORF.__init__(self, isoname, gene)


class CloneORF(SequencedORF):
    """Represents a cloned isoform. Can have assoc. PPI/PDI/func. data."""
    # self.src - source collection in hORFeome
    # self.orfid
    # self.ppis = set() # holds set of ppi objects, cases where CL_ORF is anchor of iso-matrix
    # self.ppis_as_ptr = set() # holds ppi_obj, but cases where this CL_ORF is on the partner side
    # self.has_pos_ppi = None # boolean orf has 1+ positive ppi in ppis attr.
    pass


class PacBioORF(SequencedORF):
    """Represents a polished PB transcript."""
    # self.src = 'PB'


class Exon(Biomolecule):
    def __init__(self, start, end, orf, gene):
        # self.name -> property
        # self.gene -> property
        self.orf = orf
        # self.chrom -> Biomolecule
        # self.strand -> Biomolecule
        # self.ord -> property
        self.start = start
        self.end = end
        # self.first -> property
        # self.last -> property
        # self.len -> Biomolecule
        # self.blen -> property
        # self.seq -> property
        # self.aa_seq -> property
        # self.frm_seq -> property
        # self.is_sym -> property
        self.chain = []  # list of position objects
        # self.up_exon -> property
        # self.dn_exon -> property
        self.up_junc = None
        self.dn_junc = None
        self.up_ss = None
        self.dn_ss = None
        self.cds = None  # link to cds_obj, if exists
        # self.full -> property
        Biomolecule.__init__(self)

    @property
    def name(self):
        return self.orf.name + '|E' + str(self.ord)

    @property
    def gene(self):
        return self.orf.gene

    @property
    def ord(self):
        """The ordinal of the exon, derived from position of exon in orf."""
        return self.orf.exons.index(self) + 1

    @property
    def first(self):
        """The upstream-most (i.e., first) pos_obj."""
        return self.chain[0]

    @property
    def last(self):
        """The downstream-most (i.e., last) pos_obj."""
        return self.chain[-1]

    @property
    def blen(self):
        return self.end - self.start + 1

    @property
    def seq(self):
        nt_seq = ''.join([pos.nt for pos in self.chain])
        return nt_seq

    @property
    def aa_seq(self):
        """nt-precision aa sequence. Note - three aa for codon."""
        return ''.join(pos.res.aa for pos in self.chain if pos.is_coding)

    @property
    def frm_seq(self):
        frm_chain = ''.join(str(pos.frm) for pos in self.chain)
        return frm_chain

    @property
    def is_sym(self):
        """Determine if exon is symmetrical in terms of translation frame."""
        return True if self.len % 3 == 0 else False

    @property
    def up_exon(self):
        """Get upstream exon_obj by positioning in exon list in orf_obj."""
        exon_index = self.ord - 1
        if exon_index == 0:
            # first exon in the chain
            return None
        else:
            up_exon_index = exon_index - 1
            up_exon = self.orf.exons[up_exon_index]
            return up_exon

    @property
    def dn_exon(self):
        """Get downstream exon_obj by positioning in exon list in orf_obj."""
        exon_index = self.ord - 1
        if self.ord == len(self.orf.exons):
            # last exon in the chain
            return None
        else:
            dn_exon_index = exon_index + 1
            dn_exon = self.orf.exons[dn_exon_index]
            return dn_exon

    @property
    def full(self):
        """Full representation of exon_obj."""
        ostr = ''
        ostr += self.name + ' ' + self.seq + '\n'
        ostr += self.name + ' ' + self.frm_seq + '\n'
        return ostr

    def __iter__(self):
        for pos in self.chain:
            yield pos

    def __getitem__(self, i):
        return self.chain[i]

    def get_hg38_seq(self, hg38_dict):
        """Extract the hg38-based sequences of pos in exon.  Return str.
           Note: hg38 seqs read into 0-base python dict, so decrement pos coord.
        """
        rc = {'A':'T','T':'A','C':'G','G':'C'}
        # if pos strand, return chain
        if self.strand == '+':
            hg38_chain = ''
            for pos in self.chain:
                hg38_nt = hg38_dict[self.chrom][pos.coord-1]
                hg38_chain += hg38_nt
            return hg38_chain
        # if neg strand, revcomp and return
        elif self.strand == '-':
            hg38_chain_revcomp = ''
            for nt in hg38_chain:
                nt_rc = rc[nt]
                hg38_chain_revcomp += nt_rc
            return hg38_chain_revcomp


class CDS(Exon, Biomolecule):
    """The coding (translated) portion of an Exon.

       Note - Creation is from an exon object pre-populated with pos_objs.
       The CDS object is created from a chain of residue objects and linked to
       the corresponding Exon. Assumed that incoming exon_obj has residue(s).
    """
    def __init__(self, exon, res_objs, first_pos, last_pos):
        # self.name -> property
        self.exon = exon
        # self.orf -> property
        # self.gene -> property
        # self.strand -> property
        # self.chrom -> property
        # self.ord -> property, below
        self._ord = None  # placeholder for caching (issues in resr aln looup)
        # self.start -> property
        # self.end -> property
        self.first_pos = first_pos
        self.last_pos = last_pos
        # self.first -> property, upstream-most res_obj
        # self.last -> property, downstream-most res_obj
        # self.len -> Biomolecule
        self.chain = res_objs  # includes all res, even 1/3 nt in sister exon
        # self.chain_trimmed -> property, includes major map. res only
        # self.seq -> property, from chain_trimmed
        # self.seq_all -> property, from chain
        Biomolecule.__init__(self)

    @property
    def name(self):
        return self.orf.name + '|C' + str(self.ord)

    @property
    def orf(self):
        return self.exon.orf

    @property
    def gene(self):
        return self.exon.gene

    @property
    def strand(self):
        return self.gene.strand

    @property
    def chrom(self):
        return self.gene.chrom

    @property
    def ord(self):
        """The ordinal of the cds, derived from position of cds in orf."""
        if not self._ord:
            self._ord = self.orf.cdss.index(self) + 1
        return self._ord

    @property
    def start(self):
        """Derive the start abs. coordinate, ascending on genome."""
        if self.strand == '+':
            return self.first_pos.coord
        else:
            return self.last_pos.coord

    @property
    def end(self):
        """Derive the end abs. coordinate, ascending on genome."""
        if self.strand == '+':
            return self.last_pos.coord
        else:
            return self.first_pos.coord

    @property
    def first(self):
        """Upstream-most residue in CDS."""
        return self.chain[0]

    @property
    def last(self):
        """Downstream-most residue in CDS."""
        return self.chain[-1]

    @property
    def chain_trimmed(self):
        """Chain of res_obj where this exon is major. exon (2/3, 3/3 codon)."""
        chain = []
        for res in self.chain:
            if res.exon == self.exon:  # res.exon is maj. exon, 2 or 3/3 codon
                chain.append(res)
        return chain

    @property
    def seq(self):
        """Return amino acid sequence. Only res where 2 or 3 nt in exon. This
           way, aa seq across cds of an orf can be concat. to make prot. seq.
        """
        return ''.join(res.aa for res in self.chain_trimmed)

    @property
    def seq_all(self):
        """Return amino acid sequence. Includes all res, even 1/3 mapping."""
        return ''.join(res.aa for res in self.chain)

    @property
    def full(self):
        nt_seq = ''.join([pos.nt for pos in self.exon.chain])
        ostr = '{} {}\n'.format(self.exon.name, nt_seq)
        aa_seq = ''.join([pos.res.aa for pos in self.exon.chain])
        ostr += '{} {}'.format(self.name, aa_seq)
        return ostr


class Junction(Biomolecule):
    """Represents a splice junction (exon-exon or cds-cds connection)."""
    def __init__(self, up_exon, dn_exon, orf):
        # self.name -> property
        # self.gene -> property
        # self.chrom -> Biomolecule
        # self.strand -> Biomolecule
        # self.len -> property
        self.up_exon = up_exon  # exon upstream of junction
        self.dn_exon = dn_exon  # exon downstream of junction
        self.up_ss = None  # ss_obj, repr. 5' splice site (at up_exon)
        self.dn_ss = None  # ss_obj, repr. 3' splice site (at dn_exon)
        # self.ss_seq -> property, string of splicesites (e.g. GTAG)
        # self.is_canon -> property
        # self.full -> property
        Biomolecule.__init__(self)

    @property
    def name(self):
        return self.up_exon.name + '..' + self.dn_exon.name.split('|')[-1]

    @property
    def gene(self):
        return self.orf.gene

    @property
    def len(self):
        """Define len b/c cannot use Biomolecule len, because no self.seq."""
        return len(self) # see below

    def __len__(self):
        """Length of the intron."""
        return abs(self.up_exon.last.coord - self.dn_exon.first.coord) - 1

    @property
    def ss_seq(self):
        """Splicesite dinucleotide sequences. e.g. GTAG, GCAG."""
        return self.up_ss.seq + '_' + self.dn_ss.seq

    @property
    def is_canon(self):
        """Determine if splicesites are canonical (e.g. GTAG or GCAG)."""
        # TODO: - possible bug - if dinuc. canon. but whole site not canon.
        if self.up_ss.is_canon and self.dn_ss.is_canon:
            return True
        return False

    @property
    def full(self):
        ostr = ('<' + self.up_exon.name + '>' + self.up_ss.name + '_'
                + self.dn_ss.name + '<' + self.dn_exon.name + '>')
        return ostr


class Splicesite(Biomolecule):
    """Represents the dinucleotide splice site which is a donor or accept."""
    def __init__(self, cat, ss_pos1, ss_pos2, exon, junc):
        # self.name -> property
        self.cat = cat  # donor or acceptor
        self.junc = junc
        self.exon = exon
        # self.orf -> property
        # self.gene -> property
        # self.chrom -> Biomolecule
        # self.strand -> Biomolecule
        self.ss1 = ss_pos1
        self.ss2 = ss_pos2
        # self.seq -> property
        # self.is_canon -> property
        Biomolecule.__init__(self)

    @property
    def name(self):
        return self.seq

    @property
    def orf(self):
        return self.exon.orf

    @property
    def gene(self):
        return self.exon.gene

    @property
    def seq(self):
        return self.ss1.nt + self.ss2.nt

    def is_canon(self):
        """Determine if splicesite is canonical."""
        if self.cat == 'donor' and self.seq in ['GT','GC']:
            return True
        if self.cat == 'acceptor' and self.seq in ['AG']:
            return True
        return False


class Position(Biomolecule):
    """Represents a super class of type Position (nt on the genome)."""
    def __init__(self, coord, nt, exon, orf):
        # self.name -> property
        self.exon = exon
        # self.orf -> property
        # self.gene -> property
        # self.chrom -> Biomolecule
        # self.strand -> Biomolecule
        self.coord = coord  # abs coord on chromosome
        self.nt = nt  # nucleotide (A, T, C, or G)
        # self.afrm => property
        # self.len -> Biomolecule
        # self.seq -> property
        # self.full -> property
        Biomolecule.__init__(self)

    @property
    def name(self):
        return self.seq

    @property
    def orf(self):
        return self.exon.orf

    @property
    def gene(self):
        return self.orf.gene

    @property
    def afrm(self):
        """Absolute frame, relative to the genome."""
        if self.strand == '+':
            return self.coord % 3
        elif self.strand == '-':
            afrm = self.coord % 3
            if afrm == 2:
                return 1
            if afrm == 1:
                return 2
            return 0
        else:
            return None

    @property
    def seq(self):
        return self.nt

    @property
    def full(self):
        ostr = '{} {} {}'.format(self.chrom, self.coord, self.nt)
        return ostr

    def __repr__(self):
        return self.nt


class TranscribedPosition(Position):
    """A pos_obj that is transcribed, or associated with exons of an orf. The
       pos_obj here can either be noncoding (UTR) or coding (CDS).
    """
    def __init__(self, coord, orf_idx, exon_idx, nt, frm, exon, orf, gene):
        self.idx = orf_idx  # 1-based index in orf
        # self.orf_idx -> property
        self.exon_idx = exon_idx  # 1-based index in exon
        self.frm = frm  # relative to ORF (ATG start), the frame of translation
        self.res = None  # link to residue obj
        self.is_coding = False  # turned to True if pos_obj linked to res_obj
        Position.__init__(self, coord, nt, exon, orf)

    @property
    def orf_idx(self):
        return self.idx



class SplicesitePosition(Position):
    """Represents a pos_obj that is part of the splicesite (dinucleotide)."""
    def __init__(self, coord, nt, junc, exon, orf):
        self.junc = junc
        Position.__init__(self, coord, nt, exon, orf)


class Residue(Biomolecule):
    """From three pos_objs, create a res_obj."""
    def __init__(self, orf_idx, pos_objs, exon):
        # self.name -> property
        self.codon = pos_objs  # list of three pos obj
        # self.cdss -> property, res maps to 1 or 2 cdss
        # self.cds -> property, cds res best maps (2 or 3 nt of pos)
        self.exons = self.get_exons_from_pos_objs(pos_objs)
        # self.exon  # exon res best maps (2 or 3 cod)
        # self.orf -> property
        # self.gene -> property
        # self.chrom -> Biomolecule
        # self.strand -> Biomolecule
        self.idx = self.get_aa_index(orf_idx)  # 1-based aa index, rel. to orf
        # self.nt_triplet -> property
        self.aa = self.get_translated_aa(pos_objs)
        # self.seq -> property, synonym of self.aa
        # self.p1 -> property
        # self.p2 -> property
        # self.p3 -> property
        # self.is_at_cds_edge -> property
        self.link_pos_to_this_res_obj()
        Biomolecule.__init__(self)

    @property
    def name(self):
        return str(self.idx) + '-' + self.aa

    @property
    def cdss(self):
        """All cdss that this res maps to."""
        return [exon.cds for exon in self.exons]

    @property
    def cds(self):
        """Main cds to which res maps."""
        return self.exon.cds

    @property
    def exon(self):
        """The exon to which the res maps best (3 or 2 nt of codon in exon)."""
        exons = [pos.exon for pos in self.codon]
        major_exon = Counter(exons).most_common(1)[0][0]
        return major_exon

    @property
    def orf(self):
        return self.exon.orf

    @property
    def gene(self):
        return self.orf.gene

    @property
    def nt_triplet(self):
        return ''.join(pos.nt for pos in self.codon)

    @property
    def seq(self):
        return self.aa

    @property
    def p1(self):
        return self.codon[0]

    @property
    def p2(self):
        return self.codon[1]

    @property
    def p3(self):
        return self.codon[2]

    @property
    def is_at_cds_edge(self):
        if self == self.cds.first or self == self.cds.last:
            return True
        else:
            return False

    def __repr__(self):
        return str(self.idx) + '-' + self.aa

    def full(self):
        return str(self.idx) + '-' + self.aa

    def __hash__(self):
        """Need to redefine hash func. if redefining __eq__."""
        return id(self)

    def __lt__(self, other):
        return self.idx < other.idx

    def __le__(self, other):
        return self.idx <= other.idx

    def __eq__(self, other):
        return self.idx == other.idx

    def __ne__(self, other):
        return self.idx != other.idx

    def __gt__(self, other):
        return self.idx > other.idx

    def __ge__(self, other):
        return self.idx >= other.idx

    def get_exons_from_pos_objs(self, pos_objs):
        """All exons that this res maps to, 1 or 2 exons."""
        exons = []
        for pos in self.codon:
            if pos.exon not in exons:
                exons.append(pos.exon)
        return exons

    def get_aa_index(self, nt_idx):
        """Return 1-based AA index."""
        aa_idx = int(math.ceil(float(nt_idx)/3))
        return aa_idx

    def get_translated_aa(self, pos_objs):
        return str(Seq(self.nt_triplet).translate())

    def link_pos_to_this_res_obj(self):
        # link pos to residue
        # turn on pos.is_coding flag to True
        self.p1.res = self
        self.p1.is_coding = True
        self.p2.res = self
        self.p2.is_coding = True
        self.p3.res = self
        self.p3.is_coding = True


# in isoalign, need empty res and cds objects
class EmptyResidue(Biomolecule):
    """Shell of a residue, expected to be untethered from any ORF."""
    def __init__(self, cds):
        self.name = '-'
        self.orf = EmptyORF()
        self.aa = '-'
        self.codon = [EmptyPosition(), EmptyPosition(), EmptyPosition()]
        self.idx = 0
        self.rfrm = '-'
        self.seq = ''  # return 0 (i.e., None) for length call in Biomolecule
        self.cds = cds  # expected to be an emptycds,
                        # one emptycds instance per orf pair (grp)
        self.is_at_cds_edge = None
        Biomolecule.__init__(self)

class EmptyPosition(Biomolecule):
    """Shell of a position."""
    def __init__(self):
        self.coord = 0
        self.name = '-'
        self.nt = '-'
        Biomolecule.__init__(self)

class EmptyCDS(Biomolecule):
    """Shell of a cds. Expected to be untethered."""
    def __init__(self):
        self.name = '-'
        self.ord = '-'
        self.seq = ''  # return 0 upon length call in Biomolecule
        Biomolecule.__init__(self)

class EmptyORF(Biomolecule):
    """Shell of an ORF. Expected to be unteathered. Needed for EmptyResidue
       lookup of orf.current_feat.
    """
    def __init__(self):
        self.current_feat = None
        Biomolecule.__init__(self)
