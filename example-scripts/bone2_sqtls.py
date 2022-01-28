# %%
import csv
import multiprocessing as mp
import os
import sys
from itertools import groupby, islice, product
from operator import attrgetter
from statistics import median

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biosurfer.analysis.sqtl import (
    get_event_counts, get_pblocks_related_to_junction,
    junction_has_drastic_effect_in_pair, split_transcripts_on_junction_usage)
from biosurfer.core.alignments import (Alignment,
                                       export_annotated_pblocks_to_tsv)
from biosurfer.core.constants import APPRIS, SQANTI, Strand
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (Chromosome, GencodeTranscript, Gene,
                                   Junction, PacBioTranscript, Transcript)
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LogNorm
from scipy.sparse import coo_matrix
from sqlalchemy.sql.expression import and_, or_, not_
from tqdm import tqdm

plt.switch_backend('agg')

data_dir = '../data/bone2'
output_dir = '../output/bone2'

db = Database('bone2')
session = db.get_session()

# %%
print('Loading transcript collapse mapping table...')
collapsed = dict()
with open(f'{data_dir}/hfobs_orf_refined.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        for tx in row['pb_accs'].split('|'):
            collapsed[tx] = row['base_acc']

# %%
print('Loading isoform expression table...')
expr_raw = pd.read_csv(f'{data_dir}/counts_with_cpm_norm.tsv', index_col=0, sep='\t')
expr_raw = expr_raw.drop(columns=expr_raw.columns[:13])
tx_to_idx = {tx: i for i, tx in enumerate(expr_raw.index)}
# N = len(expr_raw.index)
ijv = ((tx_to_idx[collapsed.get(tx, tx)], i, 1) for i, tx in enumerate(expr_raw.index))
row, col, data = zip(*ijv)
transform = coo_matrix((data, (row, col))).tocsr()

# %%
all_pb_ids = set(row.accession for row in session.query(PacBioTranscript.accession))
expression = pd.DataFrame(data=transform.dot(expr_raw.to_numpy()), index=expr_raw.index, columns=expr_raw.columns)
expression = expression.filter(items=all_pb_ids, axis=0).rename(lambda name: name.split('_')[0], axis=1)
for days in (0, 2, 4, 10):
    expression[f't{days}'] = expression[[col for col in expression.columns if col.startswith(f't{days}')]].mean(axis=1)
expression['total'] = expression.iloc[:, -4:].sum(axis=1)
expression['frac_total'] = expression['total'] / sum(expression['total'])

# %%
print('Loading sQTL table...')
sqtls_raw = pd.read_csv(f'{data_dir}/coloc_sig_pc_full.tsv', sep='\t', nrows=None)
sqtls = sqtls_raw[sqtls_raw['gene_type'] == 'protein_coding']
sqtls = sqtls[sqtls['gene_name'] == 'TPM2']
# sqtls = sqtls.sample(frac=0.01, axis=0)
# with open(f'{data_dir}/problems.txt') as f:
#     problems = set(line.strip() for line in f)
# sqtls = sqtls[sqtls['gene_name'].isin(problems)]
def optional_list(things):
    list_things = sorted(set(things))
    if len(list_things) == 1:
        return list_things[0]
    return list_things
sqtls = sqtls.groupby('event_id').agg(optional_list)

def stem(accession):
    return accession.split('.')[0]
def similar(acc1, acc2):
    return stem(acc1) == stem(acc2)
gene_id_stems = {stem(gene_id) for gene_id in sqtls['gene']}
gene_ids = {row.accession for row in
    session.query(Gene.accession).where(
        or_(
            False,
            *(Gene.accession.startswith(gene_id_stem) 
                for gene_id_stem in gene_id_stems)
        )
    )
}
# db.project_feature_mappings(gene_ids=gene_ids)

# %%
gene_query = session.query(Gene).join(Gene.transcripts).\
    where(
        and_(
            Transcript.sequence != None,
            Gene.accession.in_(gene_ids),
            not_(Gene.accession.contains('_', autoescape=True))
        )
    )
# print(str(gene_query))
print('Loading database objects...')
genes = {stem(gene.accession): gene for gene in gene_query}

# %%
# def abundant_and_coding(transcript: 'Transcript'):
#     return transcript.accession in over3counts and len(transcript.orfs) > 0

def sortkey(transcript: 'Transcript'):
    return -expression.total[transcript.accession], getattr(transcript, 'sqanti', SQANTI.OTHER)

def abundance(transcript: 'Transcript'):
    try:
        return expression.loc[transcript.accession]['total']
    except KeyError:
        return 0.0

def get_augmented_sqtl_record(row):
    chr, gene_id, start, end, h4 = row
    try:
        gene = genes[stem(gene_id)]
    except KeyError:
        return None
    if gene.strand is Strand.MINUS:
        start, end = end, start
    junc = Junction(start, end, chr, gene.strand)
    junc_info = {
        'chr': chr,
        'strand': str(junc.strand),
        'donor': start,
        'acceptor': end,
        'gene': gene.name,
        'H4': h4,
        'pairs': 0
    }
    pb_transcripts = {tx for tx in gene.transcripts if isinstance(tx, PacBioTranscript) and all(orf.has_stop_codon for orf in tx.orfs)}
    lacking, containing = split_transcripts_on_junction_usage(junc, filter(attrgetter('orfs'), pb_transcripts))
    containing = sorted(containing, key=sortkey)
    lacking = sorted(lacking, key=sortkey)
    if not containing or not lacking:
        return None
    
    gc_transcripts = sorted((tx for tx in gene.transcripts if isinstance(tx, GencodeTranscript)), key=attrgetter('appris'), reverse=True)
    anchor_tx = gc_transcripts[0]

    pairs = list(product(lacking, containing))
    n_pairs = len(pairs)
    alns = []
    for tx1, tx2 in pairs:
        with ExceptionLogger(f'Error for {tx1.name} and {tx2.name}'):
            alns.append(Alignment(tx1.protein, tx2.protein))
    # junc_info['pairs'] = n_pairs
    abundance_denom = sum(abundance(tx) for tx in containing) * sum(abundance(tx) for tx in lacking)
    weights = {
        (tx1, tx2): abundance(tx1) * abundance(tx2) / abundance_denom
        for tx1, tx2 in pairs
    }
    assert abs(sum(weights.values()) - 1) < 2**-8

    pblocks = get_pblocks_related_to_junction(junc, alns)

    # calculate fraction of pairs where one isoform is NMD and other is not
    junc_info['NMD'] = None
    with ExceptionLogger(f'Error for {gene.name} {junc}'):
        junc_info['NMD'] = sum(
            float(tx1.primary_orf.nmd ^ tx2.primary_orf.nmd) * weights[tx1, tx2]
            for tx1, tx2 in pairs
        )

    # calculate weighted mean of change in sequence length for junction-related pblocks
    junc_info['avg_delta_length'] = np.average(
        [pblock.delta_length for pblock in pblocks],
        weights = [weights[pblock.anchor.transcript, pblock.other.transcript] for pblock in pblocks]
    ) if pblocks else 0.0

    # classify "drastic" vs. "subtle" changes for each pair
    pair_pblocks = groupby(pblocks, key=attrgetter('anchor', 'other'))
    junc_info['drastic_effect_frequency'] = sum(
        float(junction_has_drastic_effect_in_pair(pblocks=list(pblock_group))) * weights[p1.transcript, p2.transcript]
        for (p1, p2), pblock_group in pair_pblocks
    )

    counts = get_event_counts(pblocks)
    freqs = {event.name: count/n_pairs for event, count in counts.items()}
    junc_info.update(freqs)

    if not os.path.isdir(f'{output_dir}/{gene.name}'):
        os.mkdir(f'{output_dir}/{gene.name}')

    export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks)

    force_plotting = False
    fig_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
    if force_plotting or not os.path.isfile(fig_path):
        isoplot = IsoformPlot(
            gc_transcripts + [None] + containing + lacking + [None],
            columns = {
                'novelty': lambda tx: str(getattr(tx, 'sqanti', '')),
                'GENCODE': lambda tx: tx.gencode.name if getattr(tx, 'gencode', False) else '',
                # 'CPM': abundance_str,
            }
        )
        with ExceptionLogger(f'Error plotting {gene.name} {junc}'):
            gs = GridSpec(1, 5)
            isoplot.draw_all_isoforms(subplot_spec=gs[:-1])
            isoplot.draw_frameshifts(anchor=anchor_tx)
            isoplot.draw_region(type='line', color='k', linestyle='--', linewidth=1, track=len(gc_transcripts) + len(containing) + 0.5, start=gene.start, stop=gene.stop)
            isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#f0f0f0')
            for pblock in pblocks:
                isoplot.draw_protein_block(pblock, alpha=weights[pblock.anchor.transcript, pblock.other.transcript]**0.5)
            isoplot.draw_features()
            isoplot.draw_legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
            heatmap_ax = isoplot.fig.add_subplot(gs[-1])
            all_accessions = [getattr(tx, 'accession', None) for tx in isoplot.transcripts]
            pb_accessions = [getattr(tx, 'accession', None) for tx in isoplot.transcripts if isinstance(tx, PacBioTranscript)]
            expression_sub = expression.loc[pb_accessions][['t0', 't2', 't4', 't10', 'total']]
            expression_sub = expression_sub.reindex(index=all_accessions, fill_value=0)
            expression_sub.iloc[-1, :] = expression_sub.sum(axis=0)
            heatmap = sns.heatmap(
                expression_sub,
                ax = heatmap_ax,
                xticklabels = True,
                yticklabels = False,
                linecolor = '#dedede',
                linewidths = 0.5,
                cmap = sns.cubehelix_palette(as_cmap=True, light=0.99, dark=0.0),
                cbar_kws = {'label': 'CPM'},
                norm = LogNorm(vmin=0.1, vmax=10000.0),
                mask = expression_sub.applymap(lambda x: x == 0),
                annot = expression_sub,
                fmt = '0.1f'
            )
            heatmap.set_ylabel(None)
            isoplot.fig.set_size_inches(16, 0.35*len(gene.transcripts))
            plt.savefig(fig_path, facecolor='w', transparent=False, dpi=200, bbox_inches='tight')
            tqdm.write('\tsaved '+fig_path)
        plt.close(isoplot.fig)
    return junc_info

rows = [(row.chr, row.gene, row.start, row.end, row.H4) for row in sqtls.itertuples()]
t = tqdm(
    map(get_augmented_sqtl_record, rows),
    desc = 'Analyzing sQTL junctions',
    total = sqtls.shape[0],
    file = sys.stdout,
    unit = 'junctions'
)
records = list(record for record in t if record)

sqtls_augmented = pd.DataFrame.from_records(records)
# display(sqtls_augmented)
print(f'Annotated {sqtls_augmented.shape[0]} sQTL junctions')
sqtls_augmented.to_csv(f'{output_dir}/coloc_sqtls_annotated.tsv', sep='\t', index=False)

# %%
