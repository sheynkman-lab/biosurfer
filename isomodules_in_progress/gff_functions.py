# 160703
# GFF functions

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
