
# gff utility to extract nt sequence from genome
def extract_nt_seq(ref_genome, gff):
    """
    Using exon or cds as base feature, returns fasta of hg38-extracted nt sequence.  All
    strands are forward.  hg38 canon only.  Assume that the gff is sorted ascending.
    """
    hg = open(ref_genome).readlines()
    chr = [x.split(" ")[0][1:] for x in hg[0:50:2]]
    seq = hg[1:50:2]
    d = dict(zip(chr, seq))
    d_seq = {}
    enst_strand = {}
    # Sort gff before start, ascending start index.
    orig_gff = [x.split("\t") for x in open(gff).readlines() if x[0] != "#"]
    sorted_gff = ["\t".join(x) for x in sorted(orig_gff, key = lambda x: (int(x[3]), x[0]))]
    for line in sorted_gff:
        if line[0] == "#":
            continue
        fields = line.split("\t")
        chr = fields[0]
        start = int(fields[3]) - 1
        end = int(fields[4])
        strand = fields[6]
        if fields[2] == "CDS":
            enst = line.split("Parent=")[1].split(";")[0]
            enst_strand[enst] = strand # mark negative strands for later revcomps
            if enst not in d_seq:
                d_seq[enst] = ""
            d_seq[enst] += d[chr][start:end]
    ofile = open("ensts_extracted_cds_nt.fa", "w")
    ofile2 = open("ensts_extracted_cds_translated_aa.fa", "w")
    for k, v in d_seq.items():
        if enst_strand[k] == '-':
            v = revcomp(v)
        ofile.write(">" + k + "\n" + v + "\n")
        ofile2.write(">" + k + "\n" + translate(v).strip("_") + "\n")


# holds annotation objects
# e.g. TF family, expression, GO data
class Annotation():
    """Represents a flexible class that holds annotation info."""
    pass


class Expression():
    def __init__(self, expr_dict):
        self.expr_dict = expr_dict # tissue -> rpkm
        self.avg_expr = self.compute_avg_expr()

    def __getitem__(self, tiss):
        # when expr_obj fetch by key (tissue), return value
        return self.expr_dict[tiss]


    def compute_avg_expr(self):
        tot = 0
        for k, v in self.expr_dict.items():
            v = float(v)
            tot += v
        avg_expr = tot/len(self.expr_dict)
        return avg_expr



#%%


# print out a character representation of isoforms


import os
from isomodules import isocreate
from isomodules import isofunc
import itertools
import math

data_dir = '/Users/gloriasheynkman/Documents/research_drive/projects/biosurfer/data'
odir = './results'

if not os.path.exists(odir):
	os.mkdir(odir)

# filepaths
path_gc_gtf = os.path.join(data_dir, 'gencode.v38.annotation.gtf')
path_gc_fa = os.path.join(data_dir, 'gencode.v38.pc_transcripts.fa')
path_hg38_fa = os.path.join(data_dir, 'GRCh38_canon.fa')

# load data
orf_seqs = isofunc.gc_fasta_to_orf_seq_dict(path_gc_fa)
hg38_dict = isofunc.load_hg38(path_hg38_fa)

# make gene dict for two genes
genes = ['NFYA', 'PAX5']
d_gc = isocreate.init_gen_obj_gc(path_gc_gtf, gene_list=genes)
d_gc_w_seqs = isocreate.create_and_link_seq_related_obj(d_gc, orf_seqs)
d_gc_w_seqs_w_juncs = isocreate.create_and_link_junct_and_ss_objs(d_gc_w_seqs, hg38_dict)
d_gc = d_gc_w_seqs_w_juncs

# print character iso-image file
def print_iso_char_image(orf1, orf2, scale=30):
    # helper functions
    def is_pos_in_orf(orf, i):
        for exon in orf.exons:
            for pos in exon.chain:
                if i == pos.coord:
                    return True
        return False

    def return_cat(in_orf1, in_orf2):
        if in_orf1 and in_orf2:
            return 'C'
        elif in_orf1:
            return '1'
        elif in_orf2:
            return '2'
        else:
            return '0'

    def jump_index(orf1, orf2, i):
        # set index to left-most-closest coord. position
        smallest_coord = 1000000000000
        for exon in orf1.exons + orf2.exons:
            for pos in exon.chain:
                if pos.coord > i and pos.coord < smallest_coord:
                    smallest_coord = pos.coord
        return smallest_coord

    def split_blocks(in_string):
        # split basd on contiguous char.
        split_string = [''.join(v) for k, v in itertools.groupby(in_string)]
        return split_string

    def downsample_block(size, scale):
        # downsample size of exon block by amount
        small_size = int(math.ceil(float(size)/scale))
        return small_size

    def convert_blocks_to_iso_chains(blocks, scale):
        first = ''
        second = ''
        for block in shrunk_blocks:
            for b in block:
                if b == 'C':
                    first += 'X'
                    second += 'X'
                elif b == '1':
                    first += 'X'
                    second += '-'
                elif b == '2':
                    first += '-'
                    second += 'X'
                elif b == '0':
                    first += '-'
                    second += '-'
        return first, second

    # print(out an intron-squeezed isoform representation)
    # start with lowest coord, ascend and set category
    i = 0 # current coordinate
    i_orf1 = 0 # current index of exon obj for orf1
    i_orf2 = 0 # current index of exon obj for orf2
    # first, find the lowest abs. coordinate between two orfs
    coord1 = orf1.exons[0].chain[0].coord
    coord2 = orf2.exons[0].chain[0].coord
    if coord1 == coord2:
        i = coord1
        cat = 'C'
    elif coord1 < coord2:
        i = coord1
        cat = '1'
    elif coord2 < coord1:
        i = coord2
        cat = '2'
    chain = ''
    chain += cat
    for j in range(1000):
        i += 1
        in_orf1 = is_pos_in_orf(orf1, i)
        in_orf2 = is_pos_in_orf(orf2, i)
        cat = return_cat(in_orf1, in_orf2)
        chain += cat
        if cat == '0':
            i = jump_index(orf1, orf2, i) - 1
        if i == 1000000000001:
            break
    blocks = split_blocks(chain)
    shrunk_blocks = []
    for block in blocks:
        cat = block[0]
        size = downsample_block(len(block), scale)
        shrunk_blocks.append(cat*size)
    first, second = convert_blocks_to_iso_chains(blocks, scale)
    print(first)
    print(second)

for gname, gene in d_gc.items():
    print(gname)
    print_iso_char_image(pair[0], pair[1])