# functions that wrangle (reshape, regroup, restructure) iso-objects


from isomodules import isoclass
from isomodules import isogroup



def set_gencode_gene_orfs_as_repr_orfs_from_same_prot_seq_grps(gd):
    """For all genes in the gene object dictionary, point gene to orfs that
       are representative orfs (from same-prot-seq groups). In the process,
       set gene.redundant_seq_orfs to those orfs which were filtered out
       (because it was a duplicate protein-seq but not the selected repr. orf)
    """
    for symbol, gene in gd.items():
        grps = get_same_prot_seq_grps_of_orfs_from_gene_(gene)
        repr_orfs = set(grp.repr_orf for grp in grps)  # uniq-prot-seq orfs
        gene.redundant_seq_orfs = gene.orfs - repr_orfs
        gene.orfs = repr_orfs
    return gd

def get_same_prot_seq_grps_of_orfs_from_gene_(gen_obj):
    """Given a gen_obj dict., find groups of ORFs with the same protein seq.
       and return as same-seq-grp objects.
    """
    same_seq_orfs = gen_obj.same_seq_orfs
    grps = []  # grps of same-protein-seq orf_objs
    for orfs in same_seq_orfs:
        grp = isogroup.SameSequenceGencodeGroup(orfs)
        grps.append(grp)
    return grps
