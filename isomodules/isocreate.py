#!/usr/bin/env python
# title           :isocreate.py
# description     :Functions to create web of isoform objects.
# author          :Gloria Sheynkman
# date            :May 2nd, 2019
# version         :1
# python_version  :2.6.6
# ==============================================================================

import warnings
import itertools
from collections import defaultdict

from Bio.Seq import Seq

from . import isoclass


# *****************************************************************************
# initialize gen/orf/exon objects from gencode gtf file
def init_gen_obj_gc(path_gc_gtf, gene_list=[], upload_basic_only=False,
                                    first_n_genes=0, verbose=False, iso_list=[]):
    """Given a list of genes, extract coordinate and annot. info. to create
       gen_obj and its children (orf, exon, junc, ss).
       Return a dict. (genename->gen_obj).

       Note that objects are not yet populated with pos or res objects, so no
       coding-related features (cds, res) populated.

       Input:
          path_gc_gtf - abs. or rel. path to the gtf file
          gene_list - list of gene symbols to make gen_objs of
                        If empty, make gen_obj for all genes in gtf
          upload_basic_only - make gen_obj of gencode basic annots.
          iso_list - option to only make orfs if its isoacc in this list
        Output:
           dict - dictionary of name -> gen_obj
                   Gen_obj has orf, exon, and cds objs, each with ranges.
        Note: Assumes that features ordered by gene->transcript->exon/CDS
        Note: Input only 'protein_coding' features. 'nonsense_mediated_decay'
              is also considered protein coding, but skipping for now.
    """
    genes = {} # symbol (gene name) -> gen_obj
    features_to_parse = ['gene', 'transcript', 'CDS']
    for line in open(path_gc_gtf):
        if line_is_not_valid(line): continue
        feat, chrom, strand, start, end, accs = extract_info_from_line(line)
        symbol = parse(accs, 'gene_name')
        if not line_needs_processing(gene_list, symbol, features_to_parse,
                                     feat, accs): continue
        ensg = parse(accs, 'gene_id')
        if feat == 'gene':
            gen_obj = isoclass.GencodeGene(ensg, symbol, chrom, strand,
                                           start, end)
            genes[symbol] = gen_obj
            if verbose:
                print(symbol)
        # if request to load basic transcripts, but not basic, then skip
        is_basic = is_feat_basic(accs)
        if upload_basic_only and not is_basic: continue
        # retrieve gene object, to link ORF and exon/CDS objects to it
        gen_obj = genes[symbol]
        # create and link ORF object
        if feat == 'transcript':
            tinfo = extract_transcript_relevant_info(accs)
            (isoname, enst, ensp, level, tsl, rank, enst_idx, appris,
             cds_start_nf, cds_end_nf) = tinfo
            if this_isoform_needs_to_be_made(iso_list, isoname):
                flags = [appris, level, tsl, rank, enst_idx] # for ranking
                orf_obj = isoclass.GencodeORF(start, end, enst, ensp, appris,
                                              is_basic, flags, isoname, gen_obj,
                                              cds_start_nf, cds_end_nf)
                gen_obj.orfs.add(orf_obj)
        # create and link exon object
        # in Gencode GTF, CDS *is* coding exon; here create exon_obj; later when
        # inputting sequence will the cds_obj be derived from res_objs
        elif feat == 'CDS':
            isoname, exon_ord = extract_cds_relevant_info(accs)
            if this_isoform_needs_to_be_made(iso_list, isoname):
                orf_obj = gen_obj[isoname]
                exon_obj = isoclass.Exon(start, end, orf_obj, gen_obj)
                orf_obj.exons.append(exon_obj)
    qc_check_ensure_all_orf_obj_exon_ranges_ascending(genes)
    return genes

def line_is_not_valid(line):
    if line.startswith('#') or 'ENSGR0' in line or 'ENSTR0' in line:
        return True
    else:
        return False

def extract_info_from_line(line):
    """Parse line and return data values."""
    fds = line.rstrip().split('\t')
    feat, chrom, strand, start, end = (fds[2], fds[0], fds[6],
                                       int(fds[3]), int(fds[4]))
    accs = fds[8] # accession line
    return feat, chrom, strand, start, end, accs

def line_needs_processing(gene_list, symbol, feat_list, feat, accs):
    """Determine if feature (line in gtf file) needs processing."""
    if (gene_was_requested_to_be_made(symbol, gene_list) and
           feature_needed_for_obj_creation(feat, feat_list) and
           feature_is_protein_coding(feat, accs)):
        return True
    else:
        return False

def gene_was_requested_to_be_made(gene, gene_list):
    """Determine if a gene or feature (transcript, CDS) should be created."""
    if gene_list == []:
        return True
    # Otherwise, only go forward with genes in the to-order gene list.
    if gene in gene_list:
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

def gene_not_part_of_exclusive_list(gene_list, symbol):
    """Determine if a gene_list was passed in, which lists genes for which
       to make gen_objs.  If so, determine if symbol is in the list.
       We only want to process genes in gene_list.
    """
    if len(gene_list) == 0:
        return False
    if len(gene_list) > 0 and (symbol not in gene_list):
        return True
    else:
        return False

def opted_to_load_basic_only_but_feat_not_basic(upload_basic_only, accs):
    """Determine if desire to load basic only and this feat. is not basic."""
    is_basic = is_feat_basic(accs)
    if upload_basic_only and (not is_basic):
        return True
    else:
        return False

def is_feat_basic(accs):
    """Determine if feature has basic annotation tag."""
    is_basic = True if 'tag "basic"' in accs else False
    return is_basic

def is_transcript_protein_coding(accs):
    """Determine if transcript is protein coding from transcr. feat. line."""
    tbiotype = parse(accs, 'transcript_type')
    # 'nonsense_mediated_decay' also considered prot. cod., skip for now
    # take only 'protein coding' transcripts
    if tbiotype == 'protein_coding':
        return True
    else:
        return False

def extract_transcript_relevant_info(accs):
    """Extract info from transcript line.  Return as list.

       Note: Gencode v29 had gene_status and transcript_status tags. These are
             deprecated in Gencode v30.
    """
    tags = make_acc_dict(accs) # tag -> value
    isoname = tags['transcript_name']
    enst = tags['transcript_id']
    ensp = tags['protein_id']
    level = convert_level_to_numbers(tags['level'])
    transcript_support_level = convert_tsl_to_numbers(tags)
    rank = int(isoname.split('-')[-1]) # isoform name suffix (e.g. 201)
    enst_idx = int(enst[4:].split('.')[0])
    appris = get_appris_status(accs)
    cds_start_NF = get_status_of_cds_start(tags)
    cds_end_NF = get_status_of_cds_end(tags)
    tinfo = (isoname, enst, ensp, level, transcript_support_level, rank,
             enst_idx, appris, cds_start_NF, cds_end_NF)
    return tinfo

def this_isoform_needs_to_be_made(iso_list, isoname):
    """THere is an option to only make ORFs if isoacc in list."""
    if iso_list:
        if isoname in iso_list:
            return True
        else:
            return False
    return True

def make_acc_dict(accs):
    """Convert gencode accession line into a dictionary."""
    d = defaultdict(list) # tag -> value
    pairs = accs.split(';')
    for pair in pairs:
        if len(pair) > 1:
            key, val = pair.split()
            val = val.strip('"')
            if key == 'tag': # 'tag' has multiple values, so store in list
                d['tag'].append(val)
            else:
                d[key] = val
    return d

def convert_level_to_numbers(level):
    """Level is 1, 2, 3, or NA.
       For cases where level is NA or missing, assign as 4 (lowest rank). This
       is to allow for alphanumeric sorting of isoforms
    """
    if level == 'NA':
        level = 4
    else:
        level = int(level)
    return level

def convert_tsl_to_numbers(tags):
    """Transcript_support_level (tsl) is 1, 2, 3, 4, 5, 6, NA, or missing.
       For cases where transcript_support_level is NA or missing, assign as 6
       (lowest rank). This allows for alphanumeric sorting of isoforms later
       when defining the reference transcript.
    """
    try:
        tsl = int(tags['transcript_support_level'])
    except:
        tsl = 6
    return tsl

def get_appris_status(accs):
    """Retrieve Appris status, if available.
       Note - return empty string if not appris, because later used for sorting.
    """
    appris = ''
    for acc in accs.split(';'):
        if 'appris' in acc:
            appris = acc.split('"')[1]
            appris = make_appris_names_compatible_for_alphanum_sort(appris)
            return appris
    return 'na'

def make_appris_names_compatible_for_alphanum_sort(appris):
    """Convert appris names for alphanum. sort, to select appris orf.
       Names are - appris_principal_1, appris_alternative_1, so appris_alt.
       would be listed first. Therefore, convert name to appris_secondary_*
    """
    return appris.replace('alternative', 'secondary')

def get_status_of_cds_start(tags):
    if 'cds_start_NF' in tags['tag']:
        return True
    else:
        return False

def get_status_of_cds_end(tags):
    if 'cds_end_NF' in tags['tag']:
        return True
    else:
        return False

def extract_cds_relevant_info(accs):
    """Extract info from CDS line.

       exon_ord = Exon ordinal assigned by gencode (upstream-most is 1)
    """
    isoname = parse(accs, 'transcript_name')
    exon_ord = int(accs.split('exon_number ')[1].split(';')[0]) # exon ordinal
    return isoname, exon_ord

def qc_check_ensure_all_orf_obj_exon_ranges_ascending(genes):
    """Ensure exons populated in gen_obj have proper coordinate ordering.
       Expected that intra-exon coordinates (start, end) always ascending.
       Expect that exon_obj in upstream-to-downstream ordering.
    """
    for name, gene in genes.items():
        for orf in gene.orfs:
            if orf.name == 'SHOX': return 'pass' # weirdo gencode case
            test_each_pair_of_exon_coords_ascending(orf)
            test_exon_ordering_in_orf(orf)

def test_each_pair_of_exon_coords_ascending(orf):
    """Note - found case in which the last CDS exon is 1 nucleotide (4 nt with
              stop codon)
    """
    for exon in orf.exons:
        assert exon.start <= exon.end, '{} coords not ascending'.format(exon)

def test_exon_ordering_in_orf(orf):
    # assumes that all exon ranges in an orf are non-overlapping
    starts = [exon.start for exon in orf.exons]
    if orf.strand == '+':
        test = (sorted(starts) == starts)
        assert test, '{} pos. orf, exons not read in asc.'.format(orf)
    elif orf.strand == '-':
        test = (sorted(starts, reverse=True) == starts)
        assert test, '{} neg. orf, exons not read in desc.'.format(orf)


# *****************************************************************************
# populate sequence-related objects (cds, pos, seq) into gen_obj tree
def create_and_link_seq_related_obj(gen_obj_dict, orf_seqs, verbose=False):
    """Traverse the gen_obj and it's sub-objects (orf, exon) and populate with
       sequence-related objects (pos, res, cds). Also, link them together.

       First, empty exons are populated with pos_objs, which actuates formation.
       of res_obj and, finally, cds_objs. From 3 translated pos_objs, res_obj
       is created. Then cds_objs created from blocks of res_obj that correspond
       similarly to sister exon_obj.

       Notes: If there is an overhang of nt beyond the CDS-based block lengths,
              then trim. If the trim is more than 3 nt, report error.
    """
    total_num_genes = len(gen_obj_dict)
    gene_ordinal = 0
    for gene_name, gene in gen_obj_dict.items():
        if verbose:
            gene_ordinal += 1
            if gene_ordinal % 100 == 0:
                print('pop seq: ' + gene_name + ' {} of {}'.format(gene_ordinal, total_num_genes))
        for orf in gene:
            orf_seq = orf_seqs[orf.name]
            orf_idx = 1 # 1-base orf index
            orf_seq = trim_orf_seq_if_up_to_3_nt_overhang(orf, orf_seq)
            frm_generator = itertools.cycle([0, 1, 2]) # releases pos_obj transl. frames
            orf.seq = orf_seq
            # with orf.seq set, orf.tr and orf.frac_tr avail. as property
            # traverse orf.seq and iteratively add nt as Pos obj
            # when 3 nt passed, add aa as Res obj and link to 3 Pos objs
            pos_objs = [] # holds up to three pos_obj
            stop_codon_encountered = False
            for exon in orf.exons:
                coord = get_abs_coord_of_up_most_nt(exon)
                for i in range(exon.blen):
                    nt = orf.seq[orf_idx-1]
                    frm = next(frm_generator)
                    exon_idx = i + 1 # 1-base exon index
                    # for each nt, create pos_obj and populate exon_obj
                    pos_obj = isoclass.TranscribedPosition(coord, orf_idx,
                                            exon_idx, nt, frm, exon, orf, gene)
                    exon.chain.append(pos_obj) # auto. linked to orf and gene
                    pos_objs.append(pos_obj)
                    # If 3 pos_obj created, proceed to create and link res_obj
                    if frm == 2 and not stop_codon_encountered:
                        if pos_obj_triplet_is_stop_codon(pos_objs):
                            stop_codon_encountered = True
                        else:
                            # Note - res could map to 2 exons, here last exon
                            res_obj = isoclass.Residue(orf_idx, pos_objs, exon)
                        pos_objs = []
                    orf_idx += 1
                    coord = increment_one_nt_downstream(orf.strand, coord)
            # At this point, all possible res_obj created and linked to
            # exon. Make cds_obj from linked res_objs.
            for exon in orf.exons:
                if exon_is_coding(exon):
                    res_objs = get_res_obj_chain(exon)
                    first_pos, last_pos = get_first_and_last_coding_pos_objs(
                                                                exon)
                    cds_obj = isoclass.CDS(exon, res_objs, first_pos, last_pos)
                    # note that res->cds mapping done in class as properties
                    exon.cds = cds_obj
                    cds_obj.exon = exon
    return gen_obj_dict

def trim_orf_seq_if_up_to_3_nt_overhang(orf, orf_seq):
    """Some CDS sequences in Gencode are longer than the ranges in gtf. So, trim
       sequence if up to three nt longer.
    """
    block_length = orf.blen
    seq_length = len(orf_seq)
    delta = block_length - seq_length
    if delta <= 3:
        seq_trimmed = orf_seq[0:seq_length-delta]
    else:
        warnings.warn(orf.name +
            (' len delta more than 3! ...seq len: {}, block len: {}').format(
                                                    seq_length, block_length))
    return seq_trimmed

def get_abs_coord_of_up_most_nt(exon):
    if exon.strand == '+':
        return exon.start
    else:
        return exon.end

def pos_obj_triplet_is_stop_codon(pos_objs):
    nt_triplet = ''.join(pos.nt for pos in pos_objs)
    if nt_triplet in ['TAG', 'TGA', 'TAA']: return True

def increment_one_nt_downstream(strand, coord):
    """In a strand-aware manner, increment absolute coordinate downstream."""
    if strand == '+':
        return coord + 1
    elif strand == '-':
        return coord - 1

def exon_is_coding(exon):
    """Determine if there is at least one pos_obj with a linked res_obj."""
    for pos in exon:
        if pos.res:
            return True
    return False

def get_res_obj_chain(exon):
    """Get the chain of res_objs (upstream-to-downstream) assoc. with exon_obj.
       Note - gets every res_obj, include those assoc. with only one of three
              pos_objs
    """
    res_objs = []
    for pos in exon:
        if pos.res:
            res_obj = pos.res
            if res_obj not in res_objs:
                res_objs.append(res_obj)
    return res_objs

def get_first_and_last_coding_pos_objs(exon):
    """From exon, extract first (upstream) and last (downstream) protein-coding
       (i.e., pos_obj linked to a res_obj) pos_objs.
       Assumes there was an attempt to ascribe coding objects (res_obj, cds_obj)
       to the exon before encountering this function.
    """
    coding_pos_chain = get_chain_of_coding_pos_objs(exon)
    first_pos = coding_pos_chain[0]
    last_pos = coding_pos_chain[-1]
    return first_pos, last_pos

def get_chain_of_coding_pos_objs(exon):
    """Get upstream-to-downstream chain of protein-coding pos_obj, from exon."""
    pc_chain = [pos for pos in exon if pos.res]
    return pc_chain


# *****************************************************************************
# populate and link junction and splicesite objects
def create_and_link_junct_and_ss_objs(gen_obj_dict, hg38_dict):
    """Traverse the gen_obj and all exons of each orf and create and link
       junc_obj and associated ss_obj.
    """
    gd = create_and_link_junc_obj(gen_obj_dict)
    gd = create_and_link_ss_obj(gd, hg38_dict)
    return gd

def create_and_link_junc_obj(gen_obj_dict):
    """Traverse gen_obj and create junc_obj linking up each exon. Update
       relevant bio-objections. Changes made in place.
    """
    for gene_name, gene in gen_obj_dict.items():
        for orf in gene:
            # based on first N-1 exons, create and link junc_objs
            for exon in orf.exons[0:-1]:
                junc_obj = isoclass.Junction(exon, exon.dn_exon, orf)
                exon.dn_junc = junc_obj
                exon.dn_exon.up_junc = junc_obj
    return gen_obj_dict

def create_and_link_ss_obj(gen_obj_dict, hg38_dict):
    """Traverse gene dict. and create and link splicesite objects."""
    for gene_name, gene in gen_obj_dict.items():
        for orf in gene:
            # create and link ss_objs, each ss has first (up) and second (dn) nt
            for exon in orf.exons:
                # if junc exist, create 2 ipos_obj (i=Intron), and into a ss_obj
                # had to do funky identity comparison, b/c otherwise it looks
                # for __len__ of junc obj. in catch 22
                if not exon.dn_junc is None:
                    junc = exon.dn_junc
                    coord1, coord2 = get_coord_downstream_two_ss_nt(exon)
                    nt1, nt2 = extract_hg38_nt(exon, coord1, coord2, hg38_dict)
                    ss_pos1 = isoclass.SplicesitePosition(coord1, nt1,
                                                          junc, exon, orf)
                    ss_pos2 = isoclass.SplicesitePosition(coord2, nt2,
                                                          junc, exon, orf)
                    ss_obj = isoclass.Splicesite('donor', ss_pos1, ss_pos2,
                                                 exon, junc)
                    exon.dn_junc.up_ss = ss_obj
                    exon.dn_ss = ss_obj
                if exon.up_junc:
                    junc = exon.up_junc
                    coord1, coord2 = get_coord_upstream_two_ss_nt(exon)
                    nt1, nt2 = extract_hg38_nt(exon, coord1, coord2, hg38_dict)
                    ss_pos1 = isoclass.SplicesitePosition(coord1, nt1,
                                                          junc, exon, orf)
                    ss_pos2 = isoclass.SplicesitePosition(coord2, nt2,
                                                          junc, exon, orf)
                    ss_obj = isoclass.Splicesite('acceptor', ss_pos1, ss_pos2,
                                                 exon, junc)
                    exon.up_junc.dn_ss = ss_obj
                    exon.up_ss = ss_obj
    return gen_obj_dict

def get_coord_downstream_two_ss_nt(exon):
    """Extract the absolute coordinates corresponding to nt #1 and 2 (e.g. GC)
       of a splicesite.
    """
    strand = exon.strand
    if strand == '+':
        nt1_coord = exon.last.coord + 1
        nt2_coord = exon.last.coord + 2
    elif strand == '-':
        nt1_coord = exon.last.coord - 1
        nt2_coord = exon.last.coord - 2
    return nt1_coord, nt2_coord

def extract_hg38_nt(exon, coord1, coord2, hg38_dict):
    """Extract the nt corr. to hg38."""
    chrom = exon.chrom
    if exon.strand == '+':
        nt1 = hg38_dict[chrom][coord1-1]
        nt2 = hg38_dict[chrom][coord2-1]
    elif exon.strand == '-':
        nt1 = str(Seq(hg38_dict[chrom][coord1-1]).reverse_complement())
        nt2 = str(Seq(hg38_dict[chrom][coord2-1]).reverse_complement())
    return nt1, nt2

def get_coord_upstream_two_ss_nt(exon):
    """Extract the absolute coordinates corresponding to nt #1 and 2 (e.g. AG)
       of a splicesite.
    """
    strand = exon.strand
    if strand == '+':
        nt1_coord = exon.first.coord - 2
        nt2_coord = exon.first.coord - 1
    elif strand == '-':
        nt1_coord = exon.first.coord + 2
        nt2_coord = exon.first.coord + 1
    return nt1_coord, nt2_coord


# *****************************************************************************
# init. gen/orf/exon objects from creatd gtf file (e.g. sequencing, clones)
def init_gen_obj(path_gtf, gene_list=[], version='new_acc'):
    """Given a list of genes, extract coordinate info. to create
       gen_obj and its children (orf, exon, junc, ss).
       Return a dict. (genename->gen_obj).
       This function used for isoforms not from annotation, but from
       sequencing or cloning efforts.

       Note that objects are not yet populated with pos or res objects, so no
       coding-related features (cds, res) populated.

       Input:
          path_gtf - abs. or rel. path to the gtf file, created by alignment of
                     orf sequences to genome
          gene_list - list of gene symbols to make gen_objs of
                        If empty, make gen_obj for all genes in gtf
        Output:
           dict - dictionary of name -> gen_obj
                   Gen_obj has orf and exon objects
    """
    info = {} # gene -> isoacc -> [chrom, strand, <list_of_coords>]
    for line in open(path_gtf):
        chrom, feat, strand, start, end, accs = extract_info(line)
        gene, isoacc = extract_gene_and_isoacc(accs, version)
        if not gene_needs_processing(gene, gene_list): continue
        info = update_dict(info, gene, isoacc, chrom, strand, start, end)
        # if feat == 'INS':
        #     raise UserWarning('orf has insertion:' + gene + ',' + iso_acc)
    gen_objs = {}
    for gene in info:
        # using first isoacc entry as dummy, create gen_obj
        for isoacc in info[gene]:
            chrom, strand, coords = info[gene][isoacc]
            start, end = get_start_and_end_coord(coords)
            gen_obj = isoclass.SequenceBasedGene('6K', gene, chrom, strand)
            gen_objs[gene] = gen_obj
            break
        # now, create and link orf_objs
        for isoacc in info[gene]:
            chrom, strand, coords = info[gene][isoacc]
            coords.sort()
            start, end = get_start_and_end_coord(coords)
            orf_obj = isoclass.SequencedORF('6K', start, end, isoacc, gen_obj)
            gen_obj.orfs.add(orf_obj)
            # create and link exons
            if strand == '-':
                coords.reverse()
            for i, (start, end) in enumerate(coords):
                exon_obj = isoclass.Exon(start, end, orf_obj, gen_obj)
                orf_obj.exons.append(exon_obj)
    return gen_objs

def extract_info(line):
    """Parse gtf line for info."""
    wds = line.split('\t')
    chrom, feat, strand, start, end, accs = (wds[0], wds[2], wds[6],
                                             int(wds[3]), int(wds[4]), wds[8])
    return chrom, feat, strand, start, end, accs

def extract_gene_and_isoacc(accs, version):
    gene = accs.split('gene_id "')[1].split('"')[0]
    if version == 'new_acc':
        # e.g. transcript_id "ZNF|3/4|10A02";
        isoacc = accs.split('transcript_id "')[1].split('"')[0]
    elif version == 'old_acc':
        # e.g. transcript_id "ZNF|ZNF_trinity_p03..."
        isoacc = accs.split('transcript_id "')[1].split('"')[0].split('|')[1]
    else:
        print('invalid isoacc!' + accs)
    return gene, isoacc

def gene_needs_processing(gene, gene_list):
    """Determine if this gene was requested to be processed."""
    if gene_list == []:
        return True
    else:
        if gene in gene_list:
            return True
        else:
            return False

def update_dict(info, gene, isoacc, chrom, strand, start, end):
    """Update dict. holding gene/isoacc/coords info with new line info."""
    if gene not in info:
        info[gene] = {}
    if isoacc not in info[gene]:
        info[gene][isoacc] = [chrom, strand, []]
    info[gene][isoacc][2].append([start, end])
    return info

def get_start_and_end_coord(coords):
    """Get start (left genome) and end (rigth genome) coords. of orf."""
    start = coords[0][0]
    end = coords[-1][-1]
    return start, end
