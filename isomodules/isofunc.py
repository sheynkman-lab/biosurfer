#!/usr/bin/env python
#title           :isofunc.py
#description     :Helper functions for isoform-related processing.
#author          :Gloria Sheynkman
#date            :May 7th, 2019
#version         :1
#python_version  :2.6.6
#==============================================================================

from collections import defaultdict
import pickle
import glob
import os
from operator import itemgetter
from itertools import groupby
import pandas as pd


def extract_pc_genes_from_gencode_gtf(path_gc_gtf):
    """Return a list of protein-coding genenames, from a gencode gtf."""
    features_to_parse = ['gene']
    gene_list = [] # list of GC gene symbols
    for line in open(path_gc_gtf):
        if line_is_not_valid(line): continue
        feat, chrom, strand, start, end, accs = extract_info_from_gtf_line(line)
        symbol = parse(accs, 'gene_name')
        if not line_needs_processing(gene_list, symbol, features_to_parse,
                                     feat, accs): continue
        gene_list.append(symbol)
    return gene_list

def line_is_not_valid(line):
    if line.startswith('#') or 'ENSGR0' in line or 'ENSTR0' in line:
        return True
    else:
        return False

def extract_info_from_gtf_line(line):
    """Parse line and return data values."""
    fds = line.rstrip().split('\t')
    feat, chrom, strand, start, end = (fds[2], fds[0], fds[6],
                                       int(fds[3]), int(fds[4]))
    accs = fds[8] # accession line
    return feat, chrom, strand, start, end, accs

def line_needs_processing(gene_list, symbol, feat_list, feat, accs):
    """Determine if feature (line in gtf file) needs processing."""
    if (feature_needed_for_obj_creation(feat, feat_list) and
           feature_is_protein_coding(feat, accs)):
        return True
    else:
        return False

def feature_needed_for_obj_creation(feat, feat_list):
    """Deterif feature (corr. to line in gtf file) will go towards
       instantiating object.
    """
    process_feature = True if feat in feat_list else False
    return process_feature

def feature_is_protein_coding(feat, accs):
    """Determine if feature is protein-coding from accession line."""
    if feat == 'gene':
        biotype = parse(accs, 'gene_type')
    elif feat in ['transcript', 'exon', 'CDS']:
        biotype = parse(accs, 'gene_type')
        if biotype == 'protein_coding':
            biotype = parse(accs, 'transcript_type')
        else:
            biotype = 'not protein coding!'
    else:
        raise ValueError(feat + ' is not a valid feature')
    if biotype == 'protein_coding':
        return  True
    else:
        return False

def parse(acc_line, tag):
    """Extract tag info from gc gtf line."""
    return acc_line.split(tag + ' "')[1].split('"')[0]



def gc_fasta_to_orf_seq_dict(path_fasta):
    """Convert GENCODE protein-coding transcript sequence fasta to orf sequence
       dictionary.

        Output: isoname -> transcript sequence (ATG-to-stop)
       Note - found inconsistencies in Gencode block length (CDS ranges) versus
              the sequence length. Sometimes the stop-codon is included in the
              CDS range, other times it is not. Therefore, proceed to include
              range putatively that has the stop codon. And in some cases when
              creating orf_obj, need to trimm this to fit block length.
    """
    seqs = {} # Gencode isoname -> CDS sequence
    for block in open(path_fasta).read().split('>')[1:]:
        isoname, header, seq = extract_isoname_header_and_seq(block)
        if 'ENSGR' in header: continue # skip autosomal duplicate genes
        start, end = get_cds_range(header) # CDS range, 1-base
        cds_seq = extract_cds_seq_include_stop_codon(start, end, seq)
        seqs[isoname] = cds_seq
    return seqs

def extract_isoname_header_and_seq(block):
    lines = block.strip().split('\n')
    header = lines[0]
    isoname = header.split('|')[4]
    seq = ''.join(lines[1:])
    return isoname, header, seq

def get_cds_range(header):
    """Get range of coding ORF. Range to truncate UTR+CDS+UTR to CDS."""
    for word in header.split('|'):
        if word.startswith('CDS:'):
            start, end = [int(x) for x in word[4:].split('-')]
            return start, end

def extract_cds_seq_include_stop_codon(start, end, seq):
    """Extract CDS sequence, without stop codon.
       Note - many (20K/99K) transcripts without valid stop codon. This is
              annot. under the tag 'cds_end_NF', which is added to orf_obj.
    """
    cds_seq = seq[start-1:end] # include stop codon
    return cds_seq



def load_hg38(path_hg38):
    """Load canonical chromosomes from hg38 fasta into dict."""
    hg38 = {} # chr -> seq
    for block in open(path_hg38).read().split('>')[1:]:
        lines = block.split('\n')
        chrom = lines[0].split()[0]
        if chrom.startswith('chr'):
            seq = ''.join(lines[1:])
            hg38[chrom] = seq
    return hg38



def oc_fasta_to_orf_seq_dict(path_fasta, version='new_acc'):
    """Convert fasta file to a isoname->seq dictionary."""
    orf_seqs = {} # isoname -> orf_seq
    for block in open(path_fasta).read().split('>')[1:]:
        header, sequence = block.strip().split('\n')
        if version == 'new_acc':
            isoacc = header.strip()
        elif version == 'old_acc':
            isoacc = header.rstrip().split('|')[1]
        else:
            print('invalid version number ' + path_fasta)
        orf_seqs[isoacc] = sequence
    return orf_seqs



def write_gene_pickles_to_dir(gen_dict, dname, verbose=False):
    """For all gen_obj in dict, make pickle and write-out to <dname>.
       Input:
          gen_dict - name -> gen_obj dictionary
          dname - output directory (gen_obj will be ind. pickles)
          verbose - option to print(out gen_obj being dumped, set to 'v')
    """
    make_dir_if_not_exist(os.path.join(dname))
    for name, gen_obj in gen_dict.items():
        pname = name + '.p'
        ppath = dname + pname
        if verbose:
            print('writing pickle of gen_obj: ' + ppath)
        pickle.dump(gen_obj, open(ppath, 'w'))

def make_dir_if_not_exist(dname):
    if not os.path.exists(dname):
        os.mkdir(dname)



def load_gene_pickles_to_dict(dname, verbose=False, genes=[]):
    """Load all gen_obj pickles listed in 'dname' into a dict.  Return dict.
       Input:
          dname - dir name where gene pickles are stored
          genes - if desired, select genes to upload, rather than all gen_obj
       Output:
          dict of gen_obj (name -> gen_obj)
    """
    d = {}
    dname = format_directory_string(dname)
    for ppath in glob.glob(dname + '*'):
        if verbose:
            print('loading: ' + ppath)
        pname = os.path.basename(ppath)
        name = pname.split('.')[0]
        assert name not in d, 'dup. gen_obj in pickle dir: ' + ppath
        if genes:
            if name in genes:
                d[name] = pickle.load(open(ppath))
        else:
            d[name] = pickle.load(open(ppath))
    return d



def load_clustal_alignments(path_to_alignment):
    """Load clustal pairwise alignments data into struc."""
    alns = defaultdict(list)
    for line in open(path_to_alignment).readlines()[1:]:
        gene, orf1, orf2, seq1, seq2, aln_chain = line.rstrip('\n').split('\t')
        acc = (gene, orf1, orf2)
        alns[acc] = [seq1, seq2, aln_chain]
    return alns



def load_alignment_block_calls(aln_data, path_to_alignment_blocks):
    """Load information on calls (e.g. deletion, substitution) that Sachi made on
       contiguous groups of clustal-based alignments.
       Note - ranges are 0-based
       Note - here, alignment cat. to short str (e.g. substitution to sub)
       Output: (gene, orf1, orf2) -> [[cat, start, end], etc.]
               cat examples are ins, insm, del, delm, frm, ident
               ranges (start, end) are ordered ascending
    """
    output = defaultdict(list)  # (gene, orf1, orf2) -> [cat, start, end]
    # aln_blocks: (gene, orf1, orf2) -> (start, end) -> [cat, sq1, sq2, aln]
    #                     acc        ->    range     ->        info
    aln_blocks = defaultdict(lambda: defaultdict(list))
    for line in open(path_to_alignment_blocks).readlines()[1:]:
        acc, range, info = extract_info_from_line(line)
        aln_blocks[acc][range] = info
    for acc, range_dict in aln_blocks.items():
        aln_len = get_the_alignment_length(acc, aln_data)
        blocks = get_parsed_block_info(aln_len, range_dict)
        output[acc] = blocks
    output = fillin_orfs_w_ident_seq_bc_not_repr_in_blocks(aln_data, output)
    return output

def extract_info_from_line(line):
    wds = line.rstrip('\n').split('\t')
    gene, orf1, orf2, cat, start, end, seq1, seq2, aln_chain = wds
    cat = convert_category_name(cat)  # e.g. deletion to del
    acc = (gene, orf1, orf2)
    range = (int(start), int(end))
    info = [cat, seq1, seq2, aln_chain]
    return acc, range, info

def convert_category_name(category):
    """Convert category names to a short form. For example deletion to del."""
    name_map = {'deletion': 'del', 'deletion_with_mismatch_at_edge': 'delm',
                'frameshift': 'frm', 'insertion': 'ins',
                'insertion_w_mismatch_at_edge': 'insm', 'substitution': 'sub'}
    return name_map[category]

def get_the_alignment_length(acc, aln_data):
    """For the pairwise ORF alignment, get the length of the alignment."""
    return len(aln_data[acc][0])

def get_parsed_block_info(aln_len, range_dict):
    """Input has the ranges of non-identical alignment calls.
       Supply a list of all ranges (identical blocks, non-identitcal blocks)
       so that every alignment char. is represented in the set of ranges.

       Input:
        aln_len - length of the aln, can be more than len of orf b/c of indels
        range_dict - dict. of the ranges where seq. is diff. (e.g. sub, del)

       Example: ORF of length 10
         Input: {(2,3):['deletion', 'MA', '--', 'gg']} (range_dict)
         Output: [('ID', 1, 2), ('DL', 3, 4), ('ID', 5, 10)]
    """
    ident_ranges = get_ranges_of_identical_sequence_blocks(aln_len, range_dict)
    range_dict.update(ident_ranges)
    blocks = reformat_block_calls_for_output(range_dict)  # see e.g. above
    return blocks

def get_ranges_of_identical_sequence_blocks(aln_len, range_dict):
    """Sachi only gave me ranges for which the aln. was non-identical seq.
       Therefore, heere determining ranges of samee-sequeence stretches."""
    all_idx = get_list_of_aln_idx_with_identical_sequence(aln_len, range_dict)
    ident_ranges = {}  # (st, en) -> ['ID']
    for k, g in groupby(enumerate(all_idx), lambda x: x[0] - x[1]):
        group = map(itemgetter(1), g)
        ident_ranges[(group[0], group[-1])] = ['ident']
    return ident_ranges

def get_list_of_aln_idx_with_identical_sequence(aln_len, range_dict):
    """Return a list of 0-based indices that are same-aa-sequence."""
    all_idx = set(range(aln_len))  # start with all indices
    for (start, end) in range_dict:
        sub_idx = get_list_of_aln_idx_w_nonidentical_sequence(start, end)
        all_idx = all_idx.difference(sub_idx)
    return all_idx  # at this point, all_idx are aln idx with identical seq.

def get_list_of_aln_idx_w_nonidentical_sequence(start, end):
    return set(range(start, end + 1))

def reformat_block_calls_for_output(range_dict):
    """Reformat the aln seq block calls into a list.

       Example: ORF of length 10
         Input: {(2,3):['deletion', 'MA', '--', 'gg']} (range_dict)
         Output: [('ID', 1, 2), ('DL', 3, 4), ('ID', 5, 10)]
    """
    output = []
    for (start, end) in sorted(range_dict.keys()):
        cat = range_dict[(start, end)][0]
        output.append([cat, start, end])
    return output

def fillin_orfs_w_ident_seq_bc_not_repr_in_blocks(aln_data, output):
    """Because aln block calls are only of mismatched sequences, there are no
       blocks listed for orf pairs of identical sequence. This is a problem,
       because these orf pairs are not represented in the aln_block file so
       *no* blocks get made. Therefore, here, filling in block info of full
       identical sequence.
    """
    for acc in aln_data:
        if acc not in output:
            orflen = len(aln_data[acc][0])
            output[acc] = [['ident', 0, orflen - 1]]
    return output



def load_tf_list(path_sachi, path_gloria, path_lambert):
    """Load the lists of 1300 TFs (Sachi) and 1600 TFs (Gloria)."""
    tf_sachi = [line.rstrip().split('\t')[1] for line in open(path_sachi)][1:]
    tf_gloria = [line.split('\t')[0] for line in open(path_gloria)]
    tf_lambert = []
    for line in open(path_lambert).readlines()[1:]:
        wds = line.split('\t')
        gene, is_tf = wds[1], wds[3]
        if 'Yes' in is_tf:
            tf_lambert.append(gene)
    return tf_sachi, tf_gloria, tf_lambert



def get_list_of_6k_genes(path_6k_fa):
    """Return a list of genenames in 6k collection."""
    genes = []  # list of genenames in 6k
    for block in open(path_6k_fa).read().split('>')[1:]:
        genes.append(block.split('|')[0])
    return genes



def get_updated_genenames(path_gene_map):
    """Upon upgrade to Gencode v30, several genenames were changed. Therefore
       genenames assoc. with 6K collection need to be updated to GC30.
    """
    gmap = {}  # {'orig_symbol': 'GC30_symbol'}
    for line in open(path_gene_map):
        wds = line.rstrip().split('\t')
        orig_name = wds[0]
        new_name = wds[5]
        gmap[orig_name] = new_name
    return gmap



def load_tf_genenames(fpath):
    """Load from manually corrected list of GC v30 genenames, corr. to the
       aggregate TF list (gloria + sachi + lambert).
       Notes - found that genenames as part of previously-assembled lists
               (gloria, sachi, lambert) used old GC genenames, so manually
               updated them to name listed in GCv30
    """
    genes = []
    for line in open(fpath).readlines()[1:]:
        wds = line.rstrip().split('\t')
        gene = wds[5]
        if gene == '-':
            # skip over DUX genes, which are achromosomal
            continue
        genes.append(gene)
    return genes



def load_domain_mappings(fpath_domain_mappings, fpath_pfam_names, has_eval=False):
    """Load from list of domain mappings received from Sachi.
       Also, add in the pfam domain names.

       Ouput structure:
        ensp -> [[pfam, name, cat, eval, start, end]]
        (e.g. ENSP123 -> [[pfam123, ZF, DBD, 0.0001, 4, 10]])
    """
    domains = defaultdict(list)
    pfam_names = load_pfam_names(fpath_pfam_names)
    for line in open(fpath_domain_mappings).readlines()[1:]:
        wds = line.rstrip().split('\t')
        if has_eval:
            ensp, pfam, evalue, cat, start, end = wds
        else:
            ensp, pfam, cat, start, end = wds
            evalue = -1
        try:
            domain_name = pfam_names[pfam]
        except:
            domain_name = 'unknown'
        domains[ensp].append([pfam, domain_name, cat, evalue, start, end])
    return domains

def load_pfam_names(fpath):
    names = {}  # pfam -> name
    for line in open(fpath).readlines()[1:]:
        wds = line.rstrip().split('\t')
        name, pfam = wds[0:2]
        names[pfam] = name
    return names



def load_regulatory_domain_mappings(fpath):
    """Load from a list of mappings of regulatory domain (Act, Repr, etc.)
       and the 6K isoform sequences.

       Output structure:
        isoacc -> [[acc, name, cat, start, end]]
        (e.g. ZNF123 -> [[REG00004, '-', Activator, 4, 10]])
    """
    domains = defaultdict(list)
    for line in open(fpath).readlines()[1:]:
        wds = line.rstrip().split('\t')
        isoacc, regacc, cat, start, end = wds[1], wds[2], wds[3], wds[9], wds[10]
        domains[isoacc].append([regacc, regacc, cat, start, end])
    return domains



def load_uniprot_dbd_mappings(path_domain_table, path_domain_names, non_zf_only=False):
    """Load Uniprot-domain-mapping info, thus creating a dictionary of mapped
       domains.
       Expectation is that the domain table maps each Uniprot DBD to the
       best representative isoform (first choice Gencode appris, then 6K)

       Output structure:
        gene -> isoacc -> [[acc, name, cat, start, end]]
        (e.g. ENST -> [[PFAM123, 'Myb_domain', DBD, 4, 10]])
    """
    domains = defaultdict(lambda: defaultdict(list))
    names = load_sachi_assigned_dbd_pfam_names(path_domain_names)
    for line in open(path_domain_table).readlines()[1:]:
        wds = line.rstrip().split('\t')
        cat = wds[8]
        if cat == 'DBD':
            gene, isoacc, pfam, start, end = wds[0], wds[3], wds[7], wds[12], wds[13]
            try:
                name = names[pfam]
            except:
                # pfams that are marked as DBD, but are not
                # part of the 44 sequence-specific pfams
                # confirmed that all DBDs have a name (190913)
                print(line)
            if non_zf_only:
                if 'ZF' not in name.upper():
                    domains[gene][isoacc].append([pfam, name, 'DBD', -1, start, end])
            else:
                domains[gene][isoacc].append([pfam, name, 'DBD', -1, start, end])
    return domains

def load_sachi_assigned_dbd_pfam_names(fpath):
    """Load DBD pfam names to dictionary."""
    names = defaultdict()
    for line in open(fpath).readlines()[1:]:
        wds = line.rstrip().split('\t')
        dbd_name, pfam = wds
        names[pfam] = dbd_name
    return names



def load_man_correct_gc30_tf_genenames(fpath):
    """Not all lambert/gloria/sachi 1700 tfs in GC30, so I manually corrected.
       Return a sorted list of unique TF symbols
    """
    df = pd.read_table(fpath)
    orig = set(df['orig_genename'])
    corr = set(df['new_genename'])
    tfs = orig | corr
    return tfs


def load_appris_principle_isonames(fpath):
    """Load list of appris princple isoforms for gencode tf genes."""
    return [x.rstrip().split('\t')[1] for x in open(fpath)]


def load_6k_isoacc_map(fpath):
    """Load the human readable list of 6K isoacc."""
    accs = defaultdict() # isoacc -> isoacc2 (readable), e.g. tirnity123 -> ZNF_ORF1
    for line in open(fpath).readlines()[1:]:
        isoacc, isoacc2 = line.rstrip().split('\t')[1:3]
        accs[isoacc] = isoacc2
    return accs
