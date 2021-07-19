#!/usr/bin/env python
# title           :isoalign.py
# description     :Classes that represent alignments between iso. objects.
#                  Also, functions to compute and print-out aa splice alignments.
# author          :Gloria Sheynkman
# date            :May 22nd, 2019
# version         :1
# python_version  :2.7.15
# ==============================================================================

from collections import Counter


###############################################################################
# Functions for finding aa splice alignments
def print_splice_aligned_aa_chains(pair, display='aa_default', show_rfrm=False, abacus=False):
    """Given two orf_obj, print(out genome-coord-based nt/aa alignment (i.e.)
    AA splice alignment).
        Input:
            pair - list of two orf_obj
            display - (option) 'aa_centered' or 'aa_exact'
            show_rfrm - (option) whether to show relative frame
    """
    # start of method
    # assume that orf1 (below) is longest orf
    orf1, orf2 = pair
    set_coords(pair, abacus)
    all_coords = get_list_of_coords(pair)
    orf1_coords = make_coord_dict(orf1)
    orf2_coords = make_coord_dict(orf2)
    # create orf-pair-flushed aa chain of res_obj
    # e.g. of -> chains = {orf1.name:[], orf2.name:[]} # orf-pair-aligned aa chains
    chains = make_orf_pair_aligned_aa_chain(all_coords, pair, orf1_coords, orf2_coords)
    set_rfrm_of_pos_in_orf(all_coords, orf1_coords, orf2_coords, pair)
    set_rfrm_of_res_in_chain(chains)
    if display == 'aa_default':
        ostr = return_default_display(pair, chains)
    elif display == 'aa_full':
        ostr = return_full_display(pair, all_coords, orf1_coords, orf2_coords, chains)
    ostr += '\n'
    return ostr




def return_default_display(pair, chains):
    """Return a string that displays aa-splice table."""
    orf1, orf2 = pair
    aa1 = get_aa_chain_string(chains[orf1.name])
    aa2 = get_aa_chain_string(chains[orf2.name])
    aa1f = aa1.replace('-','')  # 'flushed' aa seq
    aa2f = aa1.replace('-','')
    # aa1_frm = get_rfrm_chain_string_aa(chains[orf1.name])  # always 1s
    aa2_frm = get_rfrm_chain_string_aa(chains[orf2.name]).replace('1', ' ')
    cmp_str = derive_aa_cmp_code_string(pair, chains)
    # ranges_str = return_str_of_diff_aa_blocks_from_comparison_aa_string(aa1, aa2, cmp_str, pair, chains)
    ostr = ('{gene}\t{strand}\t{orf1}\t{aa1}\n'
            '{gene}\t{strand}\t{orf2}\t{aa2}\n'
            '{gene}\t{strand}\t{orf2}\t{aa2_frm}\n')
    ostr = ostr.format(gene=orf1.gene.name, strand=orf1.strand, orf1=orf1.name,
                aa1=aa1, orf2=orf2.name, aa2=aa2,
                aa2_frm=aa2_frm, cmp_str=cmp_str,
                aa1f=aa1f, aa2f=aa2f)
    return ostr

def get_aa_chain_string(aa_chain):
    """Return a string repr of a chain of residues (list of res_obj) that can
       contain None.  '-' represent None entries.
    """
    ostr = ''
    for res in aa_chain:
        if res:
            ostr += res.aa
        else:
            ostr += '-'
    return ostr

def get_rfrm_chain_string_aa(aa_chain):
    """Given the aa chain (spaced, b/t two orfs), get the rfrm chain."""
    ostr = ''
    for res in aa_chain:
        if res:
            ostr += res.rfrm
        else:
            ostr += '*'
    return ostr

def derive_aa_cmp_code_string(pair, chains):
    """From list of res_obj (i.e., chain), derive comparison string.
       D = deletion
       I = insertion
       S = substitution
       M = match
       Will be used for calling chunks of D/I/S/M for aa splice table.
       Assume that c1 is principle isoform.
    """
    orf1, orf2 = pair
    c1, c2 = chains[orf1.name], chains[orf2.name] # aa1_chain, aa2_chain
    ostr = ''
    for i in range(len(c1)):
        res1, res2 = c1[i], c2[i]
        if not res2:
            cat = 'D'
        elif not res1:
            cat = 'I'
        elif res1.aa == res2.aa and res1.rfrm == res2.rfrm:
            cat = 'M'
        elif res1.rfrm != res2.rfrm:
            cat = 'S'
        else:
            cat = 'U'
        ostr += cat
    return ostr


def return_str_of_diff_aa_blocks_from_comparison_aa_string(aa1_str, aa2_str, cmp_str, pair, chains):
    """From comparison string, which categorizes each aligned aa as
       a deletion, insertion, substitution, or unknown, find the ranges
       on the primary orf in which a change has occured.
       This is for the purpose of writing out a table the series of
       manipulations to go from principle to alternative protein seq.
    """
    ranges = get_ranges_of_diff_regions(cmp_str)
    blocks = get_extracted_aa_subseq_for_ranges(aa1_str, aa2_str, ranges)
    dels, ins, subs, unk = blocks
    orf1, orf2 = pair
    aa1_chain = chains[orf1.name]
    aa2_chain = chains[orf2.name]
    orf_spec_res_idx_chain1 = make_aa_chain_idx_seq(aa1_chain)
    c1 = orf_spec_res_idx_chain1
    orf_spec_res_idx_chain2 = make_aa_chain_idx_seq(aa2_chain)
    ostr = ''
    for st, en in dels:
        ostr += 'del:' + convert_idx(st,c1) + ',' + convert_idx(en,c1) + ';'
    for st, en, aa_seq in ins:
        ostr += 'ins:' + convert_ins_idx(st,c1) + ',' + aa_seq + ';'
    for st, en, aa_seq1, aa_seq2 in subs:
        ostr += 'sub:' + convert_idx(st,c1) + ',' + convert_idx(en,c1) + ',' + aa_seq1 + ',' + aa_seq2 + ';'
    for st, en, aa_seq1, aa_seq2 in unk:
        ostr += 'unk:' + convert_idx(st,c1) + ',' + convert_idx(en,c1) + ',' + aa_seq1 + ',' + aa_seq2 + ';'
    return ostr

def get_ranges_of_diff_regions(cmp_str):
    """Given a comparison string (e.g. DDDDDMMMMSSSSIIIDDDMMMM) with
       per-aa principle/alternative categorizations, return ranges
       of non-M (non-match) blocks.
       output:
        ranges - [[1,3,'D'], [4,6,'S']]
    """
    from itertools import groupby
    blocks = ["".join(grp) for num, grp in groupby(cmp_str)]
    ranges = []
    i = 1
    for block in blocks:
        cat = block[0]
        blen = len(block)
        start = i
        end = start + blen - 1
        ranges.append([start, end, cat])
        i = end + 1
    return ranges

def get_extracted_aa_subseq_for_ranges(aa1, aa2, ranges):
    """Given a brange (e.g. 1,3,S), extract the assoc. aa seqs.
       e.g. 1,3,S -> 1,3,MET,MVR
    """
    dels = []
    ins = []
    subs = []
    unk = []

    for brange in ranges:
        st, en, cat = brange
        aa1_block = aa1[st-1:en]
        aa2_block = aa2[st-1:en]
        if cat == 'D':
            dels.append([st, en])
        elif cat == 'I':
            ins.append([st, en, aa2_block])
        elif cat == 'S':
            subs.append([st, en, aa1_block, aa2_block])
        elif cat == 'U':
            unk.append([st, en, aa1_block, aa2_block])
    return dels, ins, subs, unk

def make_aa_chain_idx_seq(aa_chain):
    """Make a chain of the res.idx for an aa_chain.
       e.g. 123000456
            MAG---AGV
       For conversion of orf-aligned-aachain to orf-specific aa idx.
    """
    orf_spec_aa_idx = []
    for res in aa_chain:
        if res:
            orf_spec_aa_idx.append(res.idx)
        else:
            orf_spec_aa_idx.append(0)
    return orf_spec_aa_idx

def convert_idx(idx, orf_spec_res_idx_chain):
    """convert_orf_aligned_aa_idx_to_principle_isoform_res_idx
       Principle-Alternative-orf aligned aa chain includes inter-orf
       gaps.  Convert the idx of the orf-aligned-aa-chain to idx of
       the principle isoform
       Return idx as string for write-out.
    """
    orf_spec_idx = orf_spec_res_idx_chain[idx-1]
    return str(orf_spec_idx)

def convert_ins_idx(idx, orf_spec_res_idx_chain):
    """See function 'convert_idx'.  Convert idx of an insertion."""
    if idx == 1:
        orf_spec_idx = 0
    else:
        orf_spec_idx = orf_spec_res_idx_chain[idx-2]
    return str(orf_spec_idx)


def return_full_display(pair, all_coords, orf1_coords, orf2_coords, chains):
    """Return a string that displays aligned nt, aa, rfrm."""
    orf1, orf2 = pair
    afrm1 = get_abs_nt_frm_string(all_coords, orf1_coords)
    afrm2 = get_abs_nt_frm_string(all_coords, orf2_coords)
    ntfrm1 = get_nt_frm_string(all_coords, orf1_coords)
    ntfrm2 = get_nt_frm_string(all_coords, orf2_coords)
    nt1 = get_nt_chain_string(all_coords, orf1_coords)
    nt2 = get_nt_chain_string(all_coords, orf2_coords)
    rfrm1 = get_rfrm_chain_string(all_coords, orf1_coords)
    rfrm2 = get_rfrm_chain_string(all_coords, orf2_coords)
    aa1 = get_aa_codon_centered_chain_string(all_coords, orf1_coords)
    aa2 = get_aa_codon_centered_chain_string(all_coords, orf2_coords)

    # TODO: - complete tests
    # test_nt_feature_str_same_length_as_orf(afrm1, orf1_coords)
    test_nt_feature_str_same_length_as_orf(nt1, orf1_coords)
    test_nt_feature_str_same_length_as_orf(nt2, orf2_coords)

    ostr = '{gene}\tanchor absfrm\t{gene}|{strand}\t{orf1}\t{afrm1}\n'
    ostr += '{gene}\tanchor ntfrm \t{gene}|{strand}\t{orf1}\t{ntfrm1}\n'
    ostr += '{gene}\tanchor ntseq \t{gene}|{strand}\t{orf1}\t{nt1}\n'
    ostr += '{gene}\tanchor relfrm\t{gene}|{strand}\t{orf1}\t{rfrm1}\n'
    ostr += '{gene}\tanchor aaseq \t{gene}|{strand}\t{orf1}\t{aa1}\n'
    ostr += '{gene}\t other absfrm\t{gene}|{strand}\t{orf2}\t{afrm2}\n'
    ostr += '{gene}\t other ntfrm \t{gene}|{strand}\t{orf2}\t{ntfrm2}\n'
    ostr += '{gene}\t other ntseq \t{gene}|{strand}\t{orf2}\t{nt2}\n'
    ostr += '{gene}\t other relfrm\t{gene}|{strand}\t{orf2}\t{rfrm2}\n'
    ostr += '{gene}\t other aaseq \t{gene}|{strand}\t{orf2}\t{aa2}\n'
    ostr += '\n'
    # ostr += '{gene}\tprinc_flushd\t{gene}|{strand}\t{orf1}\t{aa}\n'
    # ostr += '{gene}\talter_flushd\t{gene}|{strand}\t{orf1}\t{aa}\n'
    ostr = ostr.format(gene=orf1.gene.name, strand=orf1.strand, orf1=orf1.name,
                aa1=aa1, orf2=orf2.name, aa2=aa2, nt1=nt1, nt2=nt2, rfrm1=rfrm1,
                rfrm2=rfrm2, ntfrm1=ntfrm1, ntfrm2=ntfrm2, afrm1=afrm1, afrm2=afrm2)
    return ostr

def get_abs_nt_frm_string(all_coords, orf_coords):
    """Given all coords, return a string of abs. nt frame."""
    ostr = ''
    for coord in all_coords:
        char = ''
        if coord in orf_coords:
            pos_obj = orf_coords[coord]
            if pos_obj.nt:
                char = str(pos_obj.afrm)
        if not char:
            char = '-'
        ostr += char
    return ostr

# TODO: - complete tests
def test_nt_feature_str_same_length_as_orf(fstr, coords):
    """Ensure that string representation of nt-related features has correct
       number of elements (matched length of orf).
    """
    for coord, pos_obj in coords.items():
        orf = pos_obj.orf
        break
    elements = fstr.replace('-', '')  # strip away fill-in characters
    value_error = 'nt feat. len mismatches orf len:{}\n{}\n{}\n'.format(orf, orf.seq, elements)
    olen = len(orf.seq)
    olen3 = olen - 3  # stop codon may be appended
    elen = len(elements)
    test = (olen == elen or olen3 == elen)
    if orf.name == 'PAX5-212':
        # case of cds_start_NF
        return
    assert test, ValueError(value_error)

def get_nt_frm_string(all_coords, orf_coords):
    """Given all coords, return a string of nt frames."""
    ostr = ''
    for coord in all_coords:
        char = ''
        if coord in orf_coords:
            pos_obj = orf_coords[coord]
            if pos_obj.nt:
                char = str(pos_obj.frm)
        if not char:
            char = '-'
        ostr += char
    return ostr

def get_nt_chain_string(all_coords, orf_coords):
    """Given all coords (inc. gaps), return a string of nt.
       '-' for gaps.
       If pos object at exon end, make lowercase.
    """
    ostr = ''
    for coord in all_coords:
        char = ''
        if coord in orf_coords:
            pos_obj = orf_coords[coord]
            if pos_obj.nt:
                if pos_obj in [pos_obj.exon.chain[0], pos_obj.exon.chain[-1]]:
                    char = pos_obj.nt.lower()
                else:
                    char = pos_obj.nt
            else:
                print('NoneType pos_obj in ' + pos_obj.orf.name)
        if not char:
            char = '-'
        ostr += char
    return ostr

def get_rfrm_chain_string(all_coords, orf_coords):
    """Return a chain of nt-resolution relative frame."""
    ostr = ''
    for coord in all_coords:
        if coord in orf_coords:
            pos_obj = orf_coords[coord]
            ostr += str(pos_obj.rfrm)
        else:
            ostr += '-'
    return ostr

def get_aa_codon_centered_chain_string(all_coords, orf_coords_dict):
    """Return a string repr. of a chain of residues, where the aa char. is
       always written directly above the 2nd (of 3) nt codons.
        Input:
            all_coords - sorted list of all coords b/t 2 orfs
            orf_coords_dict - coord_v2'd->pos_obj
    """
    ostr = ''
    for coord in all_coords:
        if coord in orf_coords_dict:
            pos_obj = orf_coords_dict[coord]
            if pos_obj == pos_obj.res.p2:
                ostr += pos_obj.res.aa
            else:
                ostr += ' '
        else:
            ostr += ' '
    return ostr

def get_frm_chain_string(all_coords, orf_coords):
    """Given all coords (inc. gaps), return a string of frm,
       '-' for gaps.
    """
    ostr = ''
    for coord in all_coords:
        if coord in orf_coords:
            pos_obj = orf_coords[coord]
            ostr += str(pos_obj.frm)
        else:
            ostr += '-'
    return ostr

def get_rfrm_chain_string_nt(all_coords, orf_coords):
    """Given all coords (inc. gaps), return a string of frm,
       '-' for gaps.
    """
    ostr = ''
    for coord in all_coords:
        if coord in orf_coords:
            pos_obj = orf_coords[coord]
            ostr += str(pos_obj.rfrm)
        else:
            ostr += '-'
    return ostr
