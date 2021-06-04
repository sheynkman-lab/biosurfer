



				
def extract_splicesite_pos_strand(genome_dict, chrom, exon1_coords, exon2_coords):
	"""Extract the left-right splicesite (dinucleotides) assoc. with junction.
	Junction is derived from pairs of exons input into this function.

	Input:
	geneome_dict - dictionary of chromosomes -> sequence (hg38 assumed)
	chrom - chromosome of input junction
	exon1_coords - pair of exon coords (e.g. [100, 150])
	exon2_coords - pair of exon coords (e.g. [200, 250])

	Output:
	splicesite - e.g. 'GTAG', 'GCAG'
	"""

	seq = genome_dict[chrom]

	# extract left splice site
	left_site = exon1_coords[1]
	left_ss = seq[left_site : left_site + 2]

	# extract right splice site
	right_site = exon2_coords[0]
	right_ss = seq[right_site - 3 : right_site - 1]

	return left_ss + right_ss


def extract_splicesite_neg_strand(genome_dict, chrom, exon1_coords, exon2_coords):
	"""Extract the left-right splicesite (dinucleotides) assoc. with junction.
	Junction is derived from pairs of exons input into this function.

	Input:
	geneome_dict - dictionary of chromosomes -> sequence (hg38 assumed)
	chrom - chromosome of input junction
	exon1_coords - pair of exon coords (e.g. [200, 250])
	exon2_coords - pair of exon coords (e.g. [100, 150])

	Note the differences in sequence extraction for negative strand.

	Output:
	splicesite - e.g. 'GTAG', 'GCAG'
	"""

	seq = genome_dict[chrom]

	# extract left splice site (upstream ss, closest to 5' end)
	left_site = exon1_coords[0]
	left_ss = seq[left_site - 3 : left_site - 1]
	left_ss = reverse_complement(left_ss)

	# extract right splice site (downstream ss, closest to 3' end)
	right_site = exon2_coords[1]
	right_ss = seq[right_site: right_site + 2]
	right_ss = reverse_complement(right_ss)

	return left_ss + right_ss



def genome_seq_to_dict(path_to_genome_fasta):
	"""Input abs. path to genome fasta, output dict. of canonical (chr1-22, X, Y) chroms.

	Input:
	path_to_genome_fasta - absolute path to the genome fasta, assume canonical chr. only

	Output:
	dictionary of chromosome sequences (chr1->seq)
	"""

	d = {}
	for block in open(path_to_genome_fasta).read().split(">")[1:]:
		lines = block.split('\n')
		chrom = lines[0].strip().split()[0]
		seq = ''.join(lines[1:])
		if chrom.startswith('chr'):
			d[chrom] = seq
	return d


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])