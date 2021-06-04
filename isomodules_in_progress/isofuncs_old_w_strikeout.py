# helper functions to manupulate data related to isoforms
import isoclass
from itertools import cycle
import subprocess
import glob
import os
import pickle
from Bio.SeqUtils import seq1, seq3 # convert from aa 3-letter-code to aa 1-letter
from pyliftover import LiftOver




# mutation-related functions

def create_mut_objs(path_clinvar):
    """Read in clinvar data to create Mutation objects.
        Input:
            path_clinvar - absolute path to clinvar variant pfile
        Output:
            set of mut_obj created from clinvar data_dir
    """
    subprocess.call('mac2unix {}'.format(path_clinvar), shell=True)
    # read in clinvar data into mut objects
    clinvar = set() # holds clinvar 'mutation' objects
    for line in open(path_clinvar).readlines()[1:]:
        fds = line.rstrip().split('\t')
        clinvar_id = fds[0]
        mut_type = fds[1]
        if mut_type != 'single nucleotide variant':
            continue
        hg = fds[16]
        if hg != 'GRCh38':
            continue
        try:
            aa_sub = fds[2].split('(p.')[1].split(')')[0]
            # if aa mapped to refseq, grab ref/alt mut
            ref_aa = seq1(aa_sub[0:3])
            alt_aa = seq1(aa_sub[-3:])
        except:
            ref_aa = None
            alt_aa = None
        symbol, cat_full, rs, medgen, dis, origin = fds[4], fds[6], fds[9], fds[12], fds[13], fds[15]
        chrom, pos, ref, alt, num_submit = fds[18], int(fds[19]), fds[21], fds[22], fds[25]
        chrom = 'chr' + chrom
        pos = int(pos)
        # input Pathogenic and Benign for now
        if 'Path' in cat_full:
            cat = 'Patho'
            mut_obj = isoclass.Mutation(ref_aa, alt_aa, clinvar_id, symbol, cat_full, cat, medgen, dis, origin, chrom, pos, ref, alt)
        elif 'Benign' in cat_full:
            cat = 'Benign'
            mut_obj = isoclass.Mutation(ref_aa, alt_aa, clinvar_id, symbol, cat_full, cat, medgen, dis, origin, chrom, pos, ref, alt)
        else:
            continue
        clinvar.add(mut_obj)
    return clinvar

def create_var_objs(path_exac):
    """Read exac data to create Variant objects.
        Input:
            path_exac - absolute path to exac data
        Output:
            set of var_obj created from exac data
    """
    print 'exac with liftover running...'
    # only do this once
    # subprocess.call('mac2unix {}'.format(path_exac), shell=True)
    # exac is based on hg19, so need to convert to hg38
    lo = LiftOver('hg19', 'hg38')
    exac = set() # holds exac 'variant' objects
    for line in open(path_exac).readlines()[1:]:
        fds = line.rstrip().split('\t')
        chrom, pos_hg19, ref, alt, snpid, filt, ac = fds[0:7]
        if filt == 'PASS' and len(ref) == 1 and len(alt) == 1:
            chrom = 'chr' + chrom
            hg38_coords = lo.convert_coordinate(chrom, int(pos_hg19))
            try:
                chrom, pos_hg38 = hg38_coords[0][0:2]
            except:
                continue
            if len(chrom) > 5: # only proceed with canon chrom
                continue
            pos = int(pos_hg38)
            ac_adj = fds[9]
            acc = chrom + '_' + str(pos) + '_' + ref + '_' + alt + '_' + snpid # make a unique exac-based acc
            flags = {'ac_adj':ac_adj, 'ac':ac, 'pos_hg19':pos_hg19}
            var_obj = isoclass.Variant('exac', chrom, pos, ref, alt, flags)
            exac.add(var_obj)
    return exac

def create_barrera_var_objs(path_bar):
    """Create variant objects from Barrera et al variants.
    Note: the only uniquely identifiable accession is the pos_acc (e.g. chr1_100001_C_T)
    """
    subprocess.call('mac2unix {}'.format(path_bar), shell=True)
    bar = set() # holds barrera variant objects
    for line in open(path_bar).readlines()[1:]:
        fds = line.rstrip().split('\t')
        src, aa_acc, gene, aa_sub, chrom, pos, ref, alt = fds
        chrom = 'chr' + chrom
        if chrom == 'chr-':
            continue
        pos_acc = chrom + '_' + pos + '_' + ref + '_' + alt
        pos = int(pos)
        flags = {'aa_acc':aa_acc, 'aa_sub':aa_sub, 'src':src}
        var_obj = isoclass.Variant('barrera', chrom, pos, ref, alt, )
        bar.add(var_obj)
    return bar
