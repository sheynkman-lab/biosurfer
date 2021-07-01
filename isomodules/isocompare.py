# functions for comparing gene/isoform structures and overlap



def get_all_matches_to_orfs(query_orf, gen_dict, match_mode='cds_coords'):
    """For a query_orf, grab all orfs to which it matches in the gen_dict.
       Input:
        query_orf - orf_obj
        gen_dict - a dict in the form of: name -> gen_obj
        match_mode - how to compare orf pairs ('cds_coords' or 'at_length')
    Output:
        matches - a list of (orf, matching_orf) pairs
    """
    matches = []
    query_genename = query_orf.gene.name
    for subject_orf in gen_dict[query_genename]:
        if query_orf_matches_the_subject_orf(query_orf, subject_orf, match_mode):
            matches.append(subject_orf)
    return matches


def query_orf_matches_the_subject_orf(query_orf, subject_orf, match_mode):
    """Helper function to determine if there is a match.  Returns orf pair."""
    if match_mode == 'cds_coords':
        if are_orfs_same_cds_coords(query_orf, subject_orf):
            return True
    elif match_mode == 'cds_intron_chain':
        if are_orfs_same_cds_intron_chain(query_orf, subject_orf):
            return True
    elif match_mode == 'at_length':
        if are_orfs_at_length_match(query_orf, subject_orf):
            return True
    else:
        raise UserWarning('unknown match_mode flag')
    return False

def are_orfs_same_cds_coords(orf1, orf2):
    """Determine if orfs have same cds chain, based on coords."""
    if len(orf1.cdss) == len(orf2.cdss):
        num_cdss = len(orf1.cdss)
        for i in range(num_cdss):
            cds1 = orf1.cdss[i]
            cds2 = orf2.cdss[i]
            if cds1.start != cds2.start or cds1.end != cds2.end:
                return False
        return True
    else:
        return False

def are_orfs_same_cds_intron_chain(orf1, orf2):
    """Determine if orfs have same internal cds chain, based on coords.
       Note that in this mode, start and end of ORF is not compared.
    """
    if len(orf1.cdss) == len(orf2.cdss):
        num_cdss = len(orf1.cdss)
        for i in range(num_cdss):
            cds1 = orf1.cdss[i]
            cds2 = orf2.cdss[i]
            if i == 0:
                if cds1.end != cds2.end:
                    return False
            elif i == (num_cdss - 1):
                if cds1.start != cds2.start:
                    return False
            else:
                if cds1.start != cds2.start or cds1.end != cds2.end:
                    return False
        return True
    else:
        return False



#TODO: - still need to test if this is correct code
def are_orfs_at_length_match(orf1, orf2):
    """Determine if orfs align from head to tail with only single nt mismatches.
    If match is 95%+ identical with alignment from head to tail, considered a
    match.
    """
    num_match = 0
    num_mismatch = 0
    if len(orf1) == len(orf2):
        for i, nt in enumerate(orf1):
            if orf1[i].nt == orf2[i].nt:
                num_match += 1
            else:
                num_mismatch += 1
        perc_ident = num_match/(float(num_match+num_mismatch))
        if perc_ident > 0.95:
            return True
    return False




# compare orf_obj against gc_sameseq_grp
def get_cds_coords_matched_gc_grps(orf, grps):
    """Return a list of gc_grps for which an orf (typically a clone-based orf) matches.
    Input:
    grps - genename -> gc_grp_obj
    Output:
    list of matched gc_grps
    """
    matched_orfs = set()
    for name, gc_grp in grps.items():
        gorf = gc_grp.repr_orf
        if are_orfs_same_cds_coords(orf, gorf):
            matched_orfs.add(gorf)
    return matched_orfs





def cluster_same_cds_coords_gc_into_gc_grps(gc_dict):
    """Given a gencode gen_obj dictionary (symbol -> gen_obj),
    returns a dictionary of symbol -> gc_grps.  orf_obj within the gen_obj
    are grouped by same-CDS coords."""
    # for each gene, put groups of same-coord orfs into a GC_Group
    #NOTE - GC_Group automatically selects 'represenative ORF' based on ranking of flags
    grps = {} # genename -> set of gc_grps
    for name, gen_obj in gc_dict.items():
        coords = make_cds_coords_centric_dict(gen_obj)
        for ocoords, orfs in coords.items():
            gc_grp = isoclass.GC_Group(orfs)
            if name not in grps:
                grps[name] = set()
            grps[name].add(gc_grp)
    return grps

def make_cds_coords_centric_dict(gen_obj):
    """Make a dictinoary where a tuple of gc_orf cds coords is the key.
    This groups all same-protein-sequence orfs into one for later ID of the repr.
    orf."""
    coords = {} # tuple of cds coords -> <list_of_orf_objs>, to find cases of 2+ orfs per same CDS coord set
    # coords - the dictionary
    # ocoords - the tuple of tuple repr. orf cds coords
    for orf in gen_obj.orfs:
        ocoords = orf.get_exon_ranges(as_tuple=True) # orf coords
        if ocoords not in coords:
            coords[ocoords] = []
        coords[ocoords].append(orf)
    return coords




def convert_grps_in_gc_grps_dict_to_repr_orf(gc_grps_dict):
    """Convert all gc_grps in gc_grps_dict to repr. orf.  Return dict.
    e.g. before: name -> set(gc_grp1, gc_grp2)
    e.g. after: name -> set(gc_orf1, gc_orf2)
    """
    converted_dict = {} # name -> set(orf_objs)
    for name, grps in gc_grps_dict.items():
        for grp in grps:
            if name not in converted_dict:
                converted_dict[name] = set()
            converted_dict[name].add(grp.repr_orf)
    return converted_dict



def compare_sequence_pair(first, second):
    """Compares nt or aa of two sequences.  Return display of mismatches/matches.
       Input can be orf or exon object, in which case the corr. seq will be compared.
       Or, input can be two sequence strings.
    """
    if isinstance(first, str) and isinstance(second, str):
        seq1 = first
        seq2 = second
    else:
        seq1 = obj1.seq
        seq2 = obj2.seq
    # TODO: find way to check that two params return ared both orf or both exon
    # seqs must be of same length
    if len(seq1) != len(seq2):
        raise UserWarning('isoclass.compare_sequence_pair method: need same-len seqs, input len {} and {}'.format(len(seq1), len(seq2)))
    match_chain = ''.join(['|' if seq1[i] == seq2[i] else '*' for i in range(len(seq1))])
    if '*' in match_chain:
        label = 'mismatches:'
    else:
        label = 'all matchd:'
    print('sequence 1:' + seq1)
    print(label + match_chain)
    print('sequence 2:' + seq2)



#TODO: - test method
def are_orfs_same_intron_chain(orf1, orf2):
    """Determine if two orf_obj have same set of intron coords.
       orf1, orf2 - two orf_objs
    """
    if orf1.get_intron_chain() == orf2.get_intron_chain():
        return True
    else:
        return False



def are_orfs_same_coords(orf1, orf2):
    """Determine if orfs have same exon chain."""
    if orf1.exons == orf2.exons:
        return True
    else:
        return False



def objs_overlap(obj1, obj2):
    """Determine if objs overlap on the genome by at least 1 nt."""
    if obj1.chrom == obj2.chrom and obj1.start <= obj2.end and obj2.start <= obj1.end:
        return True
    else:
        return False





# maybe delete these
def convert_gen_obj_dict_to_flattened_orf_obj_dict(gen_obj_dict):
    """Return a flattened gen_obj dict so it is orf.name -> orf_obj"""
    orf_obj_dict = {}
    for name, gen_obj in gen_obj_dict.items():
        for orf in gen_obj.orfs:
            orf_obj_dict[orf.name] = orf
    return orf_obj_dict
