# I updated isocreatealign and isoalign to be splice-centric alignment
# These alignment were based on Clustal and Sachi's old script.

# *****************************************************************************
# create and map alignemnts between orfs (protein-centric align, clustal)
def create_and_map_alignments(gene_dict, aln_data, aln_blocks):
    """Based on protein-centric aling, create aln_obj for res/cds/orf.
       Returns a list of grps. Gene_dict is modified in-place.

       Input:
        gene_dict - iso-obj dictionary populated with at least seq_objs
        aln_data - (gname, orf1, orf2) -> [orf1_seq, orf2_seq, aln_seq]
        aln_blocks - (gname, orf1, orf2) -> [[cat, st, en], [cat, st, en], etc.]
    !!!!!! NOTE - I turned off error print-outs for mis-matched length and
    for alnsb mapped to 2+ alnb!!!!!!!
    """
    grps = []  # group of orf pairs with populated aln_objs
    for symbol, gen_obj in gene_dict.items():
        orf_pairs = gen_obj.orf_pairs
        for pair in orf_pairs:
            info, anchor_orf, other_orf = get_full_aln_chain_info_and_anchor_other_orfs(
                                                    symbol, pair, aln_data)
            acc = (symbol, anchor_orf.name, other_orf.name)
            if orfs_not_in_gene_dict_bc_is_nmd_biotype(anchor_orf, other_orf):
                continue
            # make a group of the two orfs
            grp = isogroup.PairwiseAlignmentGroup(anchor_orf, other_orf)
            grps.append(grp)
            # instantiate the 'full' alignment (between orfs)
            alnf = isoalign.AlignmentFull(grp, anchor_orf, other_orf)
            grp.alnf = alnf
            aas1, aas2, alns = explode_and_return_alignment_tracks(info)
            if num_res_in_orf_mismatches_res_in_aln(anchor_orf, aas1,
                                                    other_orf, aas2):
                continue
            # make and link 'residue' aln_objs
            anchor_idx = 0
            other_idx = 0
            empty_cds = isoclass.EmptyCDS()  # 1 empty cds, imp. in cds retrieval
            for i in range(len(aas1)):
                aa1, aa2, aln = pop_off_first_characters(aas1, aas2, alns)
                res1, res2 = (isoclass.EmptyResidue(empty_cds),
                              isoclass.EmptyResidue(empty_cds))
                if aa1 != '-':
                    res1 = anchor_orf.res_chain[anchor_idx]
                    anchor_idx += 1
                if aa2 != '-':
                    res2 = other_orf.res_chain[other_idx]
                    other_idx += 1
                # create a residue-level alignment
                aln_obj = isoalign.AlignmentResidue(alnf, aln, res1, res2)
                alnf.chain.append(aln_obj)
            # get block info and make and link 'block' aln_objs
            blocks = aln_blocks[acc]
            for cat, start, end in blocks:
                alnb_chain = get_subset_of_alnr_based_on_range(alnf, start, end)
                alnb = isoalign.AlignmentBlock(cat, alnf, alnb_chain)
                # set downward/upward references
                alnf.blocks.append(alnb)
                for alnr in alnb.chain:
                    alnr.alnb = alnb
            # define and create 'subblock' aln_objs
            subblock_ranges = determine_subblock_ranges(alnf.chain)
            for start, end in subblock_ranges:
                alnsb_chain = get_subset_of_alnr_based_on_range(alnf, start, end)
                cds1 = get_cds_mapped_to_alnr_chain(alnsb_chain, 'anchor')
                cds2 = get_cds_mapped_to_alnr_chain(alnsb_chain, 'other')
                alnsb = isoalign.AlignmentSubblock(alnf, cds1, cds2, alnsb_chain)
                alnb = get_alnb_mapped_to_alnr_chain(alnsb_chain)
                alnsb.alnb = alnb
                alnb.subblocks.append(alnsb)
                alnf.subblocks.append(alnsb)
                for alnr in alnsb.chain:
                    alnr.alnsb = alnsb
        return grps

def get_full_aln_chain_info_and_anchor_other_orfs(symbol, pair, aln_data):
    """Retrieve the alignment track info for the pair of orfs. Note that Sachi
       randomized which one is orf1 vs orf2. Therefore check both orientations
       and return the ordering as found in her alignment file.
    """
    # two possible orderings of orf
    acc1 = (symbol, pair[0].name, pair[1].name)
    acc2 = (symbol, pair[1].name, pair[0].name)
    if acc1 in aln_data:
        anchor_orf, other_orf = pair[0], pair[1]
        info = aln_data[acc1]
    elif acc2 in aln_data:
        anchor_orf, other_orf = pair[1], pair[0]
        info = aln_data[acc2]
    else:
        print 'error {} {} {} not in aln_data'.format(symbol, pair[0], pair[1])
    return info, anchor_orf, other_orf


def orfs_not_in_gene_dict_bc_is_nmd_biotype(anchor_orf, other_orf):
    """Note - The gene_dict was instantiated with gencode orfs with strictly
              protein-codign biotype. But pairwise alignment Sachi did includes
              NMD candidates.
    """
    if anchor_orf is None or other_orf is None:
        return True
    else:
        return False


def explode_and_return_alignment_tracks(info):
    """Take string of sequence alignments and explode into list."""
    seq1, seq2, aln = info
    aas1 = list(seq1)
    aas2 = list(seq2)
    alns = list(aln)  # align. chain
    return aas1, aas2, alns


def num_res_in_orf_mismatches_res_in_aln(anchor_orf, aas1, other_orf, aas2):
    """Found cases where translated Gencode mRNA mismatches protein seq.
       See PAX5-212. Skip these cases for now.
    """
    anchor_len = len(anchor_orf.res_chain)
    anchor_aln_len = get_number_of_aas_in_aln_track(aas1)
    other_len = len(other_orf.res_chain)
    other_aln_len = get_number_of_aas_in_aln_track(aas2)
    if anchor_len != anchor_aln_len:
        pass
        # print 'Mis. len:{} {}-{}'.format(anchor_orf, anchor_len, anchor_aln_len)
        return True
    elif other_len != other_aln_len:
        pass
        # print 'Mis. len:{} {}-{}'.format(other_orf, other_len, other_aln_len)
        return True
    else:
        return False

def get_number_of_aas_in_aln_track(aln_track):
    """For example, aln_track 'AR--TE' has 4 aas."""
    return len(aln_track) - aln_track.count('-')


def pop_off_first_characters(aas1, aas2, alns):
    """Return the first character of each alignment 'track'."""
    aa1 = aas1.pop(0)
    aa2 = aas2.pop(0)
    aln_cat = alns.pop(0)
    return aa1, aa2, aln_cat


def get_subset_of_alnr_based_on_range(alnf, start, end):
    return [alnf.chain[i] for i in range(start, end+1)]


def determine_subblock_ranges(alnr_chain):
    """Determine the ranges corresponding to 'subblocks', defined as
       contiguous runs of same-anchor_cds/other_cds/aln_block_cat.
       For example:
        cds1 : 1122233333
        cds2 : 1122233444
        aln  : iiiiisssss
        i = identical, s = substitution
        -> ranges are 0-1, 2-4, 5-6, 7-9, because same cds1/cds2/aln
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


def get_cds_mapped_to_alnr_chain(alnr_chain, rank):
    """Retrieve the cds_obj the alnr are linked to.
       Based on range computed earlier, expectation is the group of alnrs
       will map to only one cds. If it maps to more than one cds, throw error.
       Note - sometimes will be mapped to no cds (because insertion/deletion)
    """
    mapped_cds = set()
    if rank == 'anchor':
        for alnr in alnr_chain:
            mapped_cds.add(alnr.res1.cds)
    elif rank == 'other':
        for alnr in alnr_chain:
            mapped_cds.add(alnr.res2.cds)
    if len(mapped_cds) > 1:
        # there are 2 distinct cdss, remove the emptycds
        mapped_cds = [cds for cds in mapped_cds if cds.name != '-']
        if len(mapped_cds) > 1:
            pass
            # print 'error - aln subblock maps to 2+ CDS ' + alnr.name
    if len(mapped_cds) == 0:
        return None
    return list(mapped_cds)[0]


def get_alnb_mapped_to_alnr_chain(alnr_chain):
    """Get the associated aln 'block' associated with the alnr of an alnsb.
       If population of alnr/alnb/alnsb was correct, expect only one alnb.
    """
    mapped_alnb = set()
    for alnr in alnr_chain:
        mapped_alnb.add(alnr.alnb)
    if len(mapped_alnb) > 1:
        pass
        # print 'error - aln subblock maps to 2+ aln block ' + str(alnr.alnf)
    if len(mapped_alnb) == 0:
        print 'no mapped alnb from alnr of an alnsb ' + str(alnr.alnf)
    return list(mapped_alnb)[0]
