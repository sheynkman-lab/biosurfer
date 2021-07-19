# helper functions for isoform-related processing

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
        domain_name = pfam_names.get(pfam, 'unknown')
        domains[ensp].append([pfam, domain_name, cat, evalue, start, end])
    return domains

def load_pfam_names(fpath):
    names = {}  # pfam -> name
    for line in open(fpath).readlines()[1:]:
        wds = line.rstrip().split('\t')
        name, pfam = wds[0:2]
        names[pfam] = name
    return names




