# %%
import multiprocessing as mp
import os
import sys
from itertools import groupby, product
from operator import attrgetter
from statistics import median

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sqlalchemy.orm import joinedload
from biosurfer.analysis.sqtl import (get_event_counts,
                                     get_pblocks_related_to_junction,
                                     junction_causes_knockdown_in_pair,
                                     split_transcripts_on_junction_usage)
from biosurfer.core.alignments import (TranscriptBasedAlignment,
                                       export_annotated_pblocks_to_tsv)
from biosurfer.core.database import db_session
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models import Chromosome, GencodeTranscript, Gene, Junction, Transcript, ORF
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from sqlalchemy.sql.expression import and_, or_
from tqdm import tqdm

# %%
data_dir = '../data/bone'
output_dir = '../output/bone'
print('Loading sQTL table...')
sqtls_raw = pd.read_csv(f'{data_dir}/coloc_sig_pc_full.tsv', sep='\t')
sqtls = sqtls_raw[sqtls_raw['gene_type'] == 'protein_coding']
def optional_list(things):
    list_things = sorted(set(things))
    if len(list_things) == 1:
        return list_things[0]
    return list_things
sqtls = sqtls.groupby('event_id').agg(optional_list)

chromosomes = Chromosome.from_names(set(sqtls['chr']))

def stem(accession):
    return accession.split('.')[0]
def similar(acc1, acc2):
    return stem(acc1) == stem(acc2)
gene_id_stems = {stem(gene_id) for gene_id in sqtls['gene']}
gene_ids = {row.accession for row in
    db_session.query(Gene.accession).where(
        or_(
            False,
            *(Gene.accession.startswith(gene_id_stem) 
                for gene_id_stem in gene_id_stems)
        )
    )
}
gene_query = db_session.query(Gene).join(Gene.transcripts).\
    where(
        and_(
            Transcript.sequence != None,
            Gene.accession.in_(gene_ids)
        )
    )
# print(str(gene_query))
print('Loading objects...')
genes = {stem(gene.accession): gene for gene in gene_query}

# %%
def get_augmented_sqtl_record(row):
    chr, gene_id, start, end, h4 = row
    try:
        gene = genes[stem(gene_id)]
    except KeyError:
        return None
    if all(isinstance(tx, GencodeTranscript) for tx in gene.transcripts):
        # skip genes for which we have no PacBio data
        return None
    junc = Junction(start, end, chromosomes[chr], gene.strand)
    junc_info = {
        'chr': chr,
        'strand': str(junc.strand),
        'donor': start,
        'acceptor': end,
        'gene': gene.name,
        'H4': h4,
        'pairs': 0
    }
    using, not_using = split_transcripts_on_junction_usage(junc, gene.transcripts)
    using = sorted(using, key=attrgetter('name'))
    not_using = sorted(not_using, key=attrgetter('name'))
    if not all((using, not_using)):
        return None
    
    pairs = list(product(not_using, using))
    n_pairs = len(pairs)
    alns = []
    for tx1, tx2 in pairs:
        with ExceptionLogger(f'Error for {tx1.name} and {tx2.name}'):
            alns.append(TranscriptBasedAlignment(tx1.protein, tx2.protein))
    junc_info['pairs'] = n_pairs

    pblocks = get_pblocks_related_to_junction(junc, alns)

    # calculate fraction of pairs where one isoform is NMD and other is not
    junc_info['NMD'] = None
    with ExceptionLogger(f'Error for {gene.name} {junc}'):
        junc_info['NMD'] = sum(tx1.primary_orf.nmd ^ tx2.primary_orf.nmd for tx1, tx2 in pairs) / n_pairs

    # calculate median change in sequence length for junction-related pblocks
    junc_info['median_delta_length'] = median(pblock.delta_length for pblock in pblocks) if pblocks else 0.0

    # classify knockdown vs. alteration for each pair
    # threshold = -not_using[0].protein.length * 2 // 5
    pair_pblocks = groupby(pblocks, key=attrgetter('anchor', 'other'))
    knockdown_pair_count = sum(junction_causes_knockdown_in_pair(pblocks=list(pblock_group)) for _, pblock_group in pair_pblocks)
    junc_info['knockdown_frequency'] = knockdown_pair_count / n_pairs

    counts = get_event_counts(pblocks)
    freqs = {event.name: count/n_pairs for event, count in counts.items()}
    junc_info.update(freqs)

    if not os.path.isdir(f'{output_dir}/{gene.name}'):
        os.mkdir(f'{output_dir}/{gene.name}')
    export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks)

    should_plot = (
        junc_info['knockdown_frequency'] > 0.9 or
        junc_info['MXIC'] >= 0.5 or
        junc_info['SIF'] >= 0.5 or
        junc_info['IR'] > 0 or
        junc_info['IX'] > 0
    )
    fig_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
    if should_plot and not os.path.isfile(fig_path):
        isoplot = IsoformPlot(using + not_using)
        with ExceptionLogger(f'Error plotting {gene.name} {junc}'):
            isoplot.draw_all_isoforms()
            isoplot.draw_frameshifts()
            isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#ffffb7')
            for pblock in pblocks:
                isoplot.draw_protein_block(pblock)
            isoplot.fig.set_size_inches(9, 0.5*len(gene.transcripts))
            plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
            tqdm.write('\tsaved '+fig_path)
        plt.close(isoplot.fig)
    return junc_info

rows = [(row.chr, row.gene, row.start, row.end, row.H4) for row in sqtls.itertuples()]
t = tqdm(
    map(get_augmented_sqtl_record, rows),
    desc = 'Analyzing sQTL junctions',
    total = sqtls.shape[0],
    file = sys.stdout
)
records = list(record for record in t if record)

sqtls_augmented = pd.DataFrame.from_records(records)
display(sqtls_augmented)
sqtls_augmented.to_csv(f'{output_dir}/coloc_sqtls_annotated.tsv', sep='\t', index=False)

# %%
