from portion.interval import Interval, Atomic, inf, closed
from portion.const import Bound
from itertools import accumulate
from typing import Collection, Tuple
import pickle
import os
from Bio import SeqIO


class IntegerInterval(Interval):
    """Represents an interval over the integers."""

    # convert all atomic intervals to open intervals, merge, then convert back to closed
    # from https://github.com/AlexandreDecan/portion/issues/24#issuecomment-604456362
    @staticmethod
    def expand(s):
        lower, upper = s.lower, s.upper
        if s.left is Bound.CLOSED:
            lower -= 1
        if s.right is Bound.CLOSED:
            upper += 1
        return Atomic(Bound.OPEN, lower, upper, Bound.OPEN)

    @staticmethod
    def reduce(s):
        if s.lower == inf or s.upper == -inf:
            return Atomic(Bound.OPEN, inf, -inf, Bound.OPEN)
        lower, upper = s.lower, s.upper
        if s.left is Bound.OPEN:
            lower += 1  
        if s.right is Bound.OPEN:
            upper -= 1
        return Atomic(Bound.CLOSED, lower, upper, Bound.CLOSED)

    def __init__(self, *intervals):
        super().__init__(*intervals)
        self._intervals = self.apply(IntegerInterval.expand).apply(IntegerInterval.reduce)._intervals
    
    def to_tuples(self):
        if self.empty:
            return ()
        return tuple((a.lower, a.upper) for a in self)

    @staticmethod
    def from_tuples(*intervals):
        return IntegerInterval(*[closed(a, b) for a, b in intervals])
    
    # Return number of integers contained in interval.
    @property
    def size(self):
        if self.empty:
            return 0
        return sum(s.upper - s.lower + 1 for s in self._intervals)


def get_start_end_from_lengths(lengths: Collection[int]) -> Tuple[Tuple[int, int]]:
    """Convert a series of interval lengths to a series of interval (start, end) coordinates.
    Coordinates are 1-indexed and inclusive."""
    starts = accumulate(lengths[:-1], initial=1)
    ends = accumulate(lengths)
    return tuple(zip(starts, ends))


# originally from https://gist.github.com/ihincks/6a420b599f43fcd7dbd79d56798c4e5a
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.
    
    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = np.array(colorsys.rgb_to_hls(*mc.to_rgb(c)))
    return colorsys.hls_to_rgb(c[0],1-amount * (1-c[1]),c[2])


def write_gene_pickles_to_dir(gen_dict, dname, verbose=False):
    """For all gen_obj in dict, make pickle and write-out to <dname>.
       Input:
          gen_dict - name -> gen_obj dictionary
          dname - output directory (gen_obj will be ind. pickles)
          verbose - option to print(out gen_obj being dumped, set to 'v')
    """
    if not os.path.exists(dname):
        os.mkdir(dname)
    for name, gen_obj in gen_dict.items():
        pname = name + '.p'
        ppath = dname + pname
        if verbose:
            print('writing pickle of gen_obj: ' + ppath)
        pickle.dump(gen_obj, open(ppath, 'w'))


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
    for rec in SeqIO.parse(path_fasta, 'fasta'):
        if 'ENSGR' in rec.id: continue # skip autosomal duplicate genes
        isoname = rec.id.split('|')[4]
        # get cds range (1-base)
        for word in rec.id.split('|'):
            if word.startswith('CDS:'):
                start, end = [int(x) for x in word[4:].split('-')]
        # extract cds sequence, with stop codon
        cds_seq = str(rec.seq[start-1:end]) # include stop codon
        seqs[isoname] = cds_seq
    return seqs


def load_domain_mappings(fpath_domain_mappings, fpath_pfam_names, has_eval=False):
    """Load from list of domain mappings received from Sachi.
       Also, add in the pfam domain names.

       Ouput structure:
        ensp -> [[pfam, name, cat, eval, start, end]]
        (e.g. ENSP123 -> [[pfam123, ZF, DBD, 0.0001, 4, 10]])
    """
    domains = defaultdict(list)
    pfam_names = {}
    for line in open(fpath_pfam_names).readlines()[1:]:
        wds = line.rstrip().split('\t')
        name, pfam = wds[0:2]
        pfam_names[pfam] = name
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