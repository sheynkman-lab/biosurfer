# classes that represent abstractions of isoform-related objects

from Bio.Seq import Seq

class Gene():
    def __init__(self, name, chrom, strand, start=0, end=0, ensg=None):
        self.muts = set()
class ORF():
    def __init__(self, cat, plot_ord, name, gene, start, end):
        self.muts = set()
        self.mgrps = set() # lists all orf-mapped mutation groups

    def derive_and_return_rel_exon_ranges(self):
        """Return a list of lists that show 1-base relative exon ranges.
        For easier viewing of coords."""
        ranges = self.get_exon_ranges()
        # get flat list of coords
        ranges_flat = [coord for pair in ranges for coord in pair]
        lowest_coord = min(ranges_flat)
        rel_ranges = []
        for pair in ranges:
            rel_pair = [coord - lowest_coord + 1 for coord in pair]
            rel_ranges.append(rel_pair)
        return rel_ranges

    def get_intron_chain(self):
        """Return list of pairs of intron coords (up to dn direction)."""
        intron_chain = []
        coords = [] # flat list of coords
        if self.strand == '+':
            for exon in self.exons:
                coords.append(exon.start)
                coords.append(exon.end)
        elif self.strand == '-':
            pairs = []
            for exon in self.exons:
                pairs.append([exon.start, exon.end])
            for pair in sorted(pairs):
                coords.append(pair[0])
                coords.append(pair[1])
        if len(coords) == 2: # single exon transcript lacks intron-chain
            return False
        elif len(coords) < 2:
            raise UserWarning('encountered orf coords less than 2 coords:'+exon.name)
        else:
            coords_trimmed = coords[1:-1]
        for i in range(0, len(coords_trimmed), 2):
            intron_pair = [coords_trimmed[i], coords_trimmed[i+1]]
            intron_chain.append(intron_pair)
        return intron_chain


class Exon():
    def __init__(self, ordinal, start, end, orf, gene):
        self.muts = set()


class Pos():
    def __init__(self, coord, idx, ridx, nt, frm, exon, orf, gene):
        self.muts = set() # list of mut objs. mapping to Pos, used to check if muts exist
        self.mgrps = set() # set of mgrps, each mgrp being composite of pos, mpos, res, mres

    def add_mut_and_spawn_and_link_mgrp(self, mut_obj):
        self.muts.add(mut_obj)
        self.link_mut_to_iso_objs(mut_obj)
        # only create mgrp if pos has exisitng res (i.e. codon is translated)
        if self.res != None:
            mpos_obj = Mutation_Pos(self, mut_obj) # link mres after init below
            mres_obj = Mutation_Residue(self, mpos_obj, mut_obj)
            mpos_obj.mres = mres_obj
            mgrp_obj = Mutation_Group(mut_obj, self, mpos_obj, self.res, mres_obj)
            self.mgrps.add(mgrp_obj)
            self.res.mgrps.add(mgrp_obj)
            self.orf.mgrps.add(mgrp_obj)

    def link_mut_to_iso_objs(self, mut_obj):
        self.exon.muts.add(mut_obj)
        self.orf.muts.add(mut_obj)
        self.gene.muts.add(mut_obj)

    def does_mut_match(self, mut_obj):
        if (self.gene.chrom == mut_obj.chrom and
            self.coord == mut_obj.pos):
            return True
        return False


class Mutation_Group():
    """Represents a composite of all objects corresponding to a mapped mutation.
    A mut_obj represents a mutatin annotation.  pos_obj across isoforms point to it.
    However, one pos_obj can point to 2+ mut_obj.  Therefore, need to create objects
    for each unique mapping.
    Note: creation of a mgrp_obj involves creation of mpos and mres objects.
    Note: in cases in which a mut_obj maps to a pos_obj lacking a res_obj (e.g. res_obj.aa = None)
    then mres of mgrp will be None as well."""

    def __init__(self, mut_obj, pos_obj, mpos_obj, res_obj, mres_obj):
        self.mut = mut_obj
        self.pos = pos_obj
        self.mpos = mpos_obj
        self.res = res_obj
        self.mres = mres_obj
        if self.res == None:
            self.has_res = False
        else:
            self.has_res = True
        # determine if clinvar-given ref nt matches orf-seq
        self.is_concordant = self.determine_if_mut_ref_and_orf_ref_concordant()

    def __repr__(self):
        gene = self.pos.gene
        orf = self.pos.orf
        pos = self.pos
        mpos = self.mpos
        mres = self.mres
        mut = self.mut
        odata = [orf.name, gene.chrom, gene.strand, pos, mpos.ref_nt, mpos.alt_nt, mres]
        # odata = [orf.name]
        # odata = [orf.name, gene.chrom, pos.coord, gene.strand, mpos.orf_ref_nt, mpos.ref_nt, self.is_concordant, mpos.alt_nt, mres, pos.idx]
        ostr = '\t'.join(map(str, odata))
        return ostr

    def determine_if_mut_ref_and_orf_ref_concordant(self):
        """Compare the mutant-given ref_nt with the nt in corresponding orf."""
        mut_ref = self.mut.ref_nt
        orf_ref = self.pos.nt
        strand = self.pos.gene.strand
        # revcomp if in negative strand
        if strand == '-':
            mut_ref = str(Seq(mut_ref).reverse_complement())
        if mut_ref == orf_ref:
            return True
        else:
            return False


#TODO: - decide if I want to encode ref/alt_nt as revcomp for neg strands
class Mutation_Pos():
    """Represents a unique Pos-to-ORF mapping."""
    def __init__(self, pos_obj, mut_obj):
        if pos_obj.gene.strand == '+':
            self.ref_nt = mut_obj.ref_nt
        else:
            self.ref_nt = str(Seq(mut_obj.ref_nt).reverse_complement())
        self.orf_ref_nt = pos_obj.nt
        # check if hg38 nt and orf-seq nt is concordant, sometimes orf seq mismatched
        if self.ref_nt == self.orf_ref_nt:
            self.is_ref_concordant = True
        else:
            self.is_ref_concordant = False
        if pos_obj.gene.strand == '+':
            self.alt_nt = mut_obj.alt_nt
        else:
            self.alt_nt = str(Seq(mut_obj.alt_nt).reverse_complement())
        self.nt = self.alt_nt
        self.pos = pos_obj
        self.mut = mut_obj
        self.res = pos_obj.res # may be None if res not translated
        self.mres = None # will link ref later

    def __repr__(self):
        ostr = 'mpos:' + self.nt
        return ostr


class Residue():
    def __init__(self, nt_idx, nt_triplet, exon, orf, gene):
        self.mgrps = set()

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


class Mutation_Residue(Residue):
    # inheriting from residue to use translate_codon method
    def __init__(self, pos_obj, mpos_obj, mut_obj):
        # self.ref_matches = None # does mapped ORF position nt match mutation-given ref
        self.pos = pos_obj
        self.mpos = mpos_obj
        self.res = pos_obj.res
        self.mut = mut_obj
        self.codon_position = pos_obj.frm # is mut in position 0, 1, or 2 of codon
        self.codon = self.set_mutation_codon(mpos_obj, pos_obj.res, pos_obj.frm)
        self.nt_triplet = self.get_nt_triplet_from_codon(self.codon) # have nt_triplet as attr for convenience
        # only return mutant aa if res is translated to aa (some res not translated b/c orf has early stop)
        if self.res.aa == None:
            self.aa = None
        else:
            self.aa = self.translate_codon(self.nt_triplet)

    def __repr__(self):
        """Format is AaG->AcG, K->T, to show codon change."""
        i = self.codon_position
        mut_trip = list(self.nt_triplet)
        mut_trip[i] = mut_trip[i].lower()
        mut_trip = ''.join(mut_trip)
        trip = list(self.res.nt_triplet)
        trip[i] = trip[i].lower()
        trip = ''.join(trip)
        ostr = trip + '-' + mut_trip + '\t' + self.res.aa + str(self.res.idx) + self.aa
        return ostr

    def set_mutation_codon(self, mpos_obj, res_obj, ordinal):
        """Determine which codon position the mutated nt is placed.
        Make a new codon holding 2 pos_obj and 1 mpos_obj.
        res_obj - normal residue object
        ordinal - which codon position the mpos_obj resides (0, 1, or 2)"""
        codon = res_obj.codon
        #TODO: ensure i'm not writing over original codon
        if ordinal == 0:
            mcodon = [mpos_obj] + codon[1:]
        elif ordinal == 1:
            mcodon = [codon[0]] + [mpos_obj] + [codon[2]]
        elif ordinal == 2:
            mcodon = codon[0:2] + [mpos_obj]
        else:
            raise UserWarning('could not assign mcodon for:' + res_obj.orf.name)
        return mcodon


class Variant():
    def __init__(self, db, chrom, pos, ref_nt, alt_nt, flags={}, name=''):
        self.db = db # e.g. clinvar, exac, hgmd, barrera
        self.chrom = chrom
        self.pos = pos
        self.ref_nt = ref_nt
        self.alt_nt = alt_nt
        self.flags = flags # optional flags
        self.name = name

    @property
    def pos_acc(self):
        return '{}_{}_{}_{}'.format(chrom, pos, ref_nt, alt_nt)

    def __repr__(self):
        ostr = self.pos_acc
        return ostr


class Mutation(Variant):
    def __init__(self, ref_aa, alt_aa, clinvar_id, symbol, cat_full, cat, medgen, dis, origin, *args):
        self.ref_aa = ref_aa # ref_aa and alt_aa exist if mut-refseq mapping exists
        self.alt_aa = alt_aa
        self.clinvar_id = clinvar_id # is also Variant.name
        self.symbol = symbol # mapped gene symbol, if exists
        self.cat_full = cat_full # benign/patho/vus, full desc.
        self.cat = cat # benign or patho
        self.medgen = medgen # medgen accessions (incl OMIM acc)
        self.dis = dis # disease discription
        self.origin = origin # germline or somatic
        Variant.__init__(self, 'clinvar', *args)

    def __repr__(self):
        ostr = ('clinvar_id-' + self.clinvar_id + ' ' + self.chrom + ':' +
                str(self.pos) + ' ' + self.ref_nt + '->' + self.alt_nt)
        return ostr
