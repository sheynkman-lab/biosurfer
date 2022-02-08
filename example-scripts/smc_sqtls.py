# %%
import csv
import multiprocessing as mp
import os
import re
import sys
from itertools import groupby, product
from operator import attrgetter, itemgetter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biosurfer.analysis.sqtl import (get_event_counts,
                                     get_pblocks_related_to_junction,
                                     junction_has_drastic_effect_in_pair,
                                     split_transcripts_on_junction_usage)
from biosurfer.core.alignments import (Alignment,
                                       export_annotated_pblocks_to_tsv)
from biosurfer.core.constants import SQANTI, AnnotationFlag, Strand
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (GencodeTranscript,
                                                Gene, Junction,
                                                PacBioTranscript, Transcript)
from biosurfer.plots.plotting import IsoformPlot, mpatches
from IPython.display import display
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from scipy.sparse import coo_matrix
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import and_, not_, or_, select
from tqdm import tqdm

plt.switch_backend('agg')

data_dir = '../data/smc'
output_dir = '../output/smc'

db = Database('smc')

# %%
print('Loading transcript collapse mapping table...')
collapsed = dict()
with open(f'{data_dir}/SMC_orf_refined.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        for tx in row['pb_accs'].split('|'):
            collapsed[tx] = row['base_acc']

# %%
print('Loading isoform expression table...')
expr_raw = pd.read_csv(f'{data_dir}/SMC_classification_TPM.txt', sep='\t', index_col=0, usecols=['isoform'] + [f'FL_TPM.S{i:02d}_3p..NEB_5p' for i in range(1, 13)])
tx_to_idx = {tx: i for i, tx in enumerate(expr_raw.index)}
# N = len(expr_raw.index)
ijv = ((tx_to_idx[collapsed.get(tx, tx)], i, 1) for i, tx in enumerate(expr_raw.index))
row, col, data = zip(*ijv)
transform = coo_matrix((data, (row, col))).tocsr()

# %%
pattern = re.compile(r'S\d{2}')
with db.get_session() as session:
    all_pb_ids = set(acc for acc in session.execute(select(PacBioTranscript.accession)).scalars())
expression = pd.DataFrame(data=transform.dot(expr_raw.to_numpy()), index=expr_raw.index, columns=expr_raw.columns)
expression = (expression.
    filter(items=all_pb_ids, axis=0).
    rename(columns=lambda name: match[0] if (match := pattern.search(name)) else name)
)
expression['FBS'] = expression.iloc[:, 0:12:2].mean(axis=1)
expression['noFBS'] = expression.iloc[:, 1:12:2].mean(axis=1)
expression['average'] = expression.iloc[:, 0:12].mean(axis=1)

# %%
def first(x):
    f = sorted(filter(None, x))
    return f[0] if f else None

print('Loading sQTL table...')
sqtls_raw: pd.DataFrame = pd.read_csv(f'{data_dir}/smc-sqtls.tsv', sep='\t', nrows=None)
sqtls_raw['gene_id'] = sqtls_raw['Cluster and ENSEMBL Gene ID'].apply(lambda x: x.split(':')[-1])
sqtls_raw['start'] = sqtls_raw['Cluster and ENSEMBL Gene ID'].apply(lambda x: int(x.split(':')[1]))
sqtls_raw['end'] = sqtls_raw['Cluster and ENSEMBL Gene ID'].apply(lambda x: int(x.split(':')[2]))
sqtls_raw['strand'] = sqtls_raw['Cluster and ENSEMBL Gene ID'].apply(lambda x: Strand.PLUS if '+' in x else Strand.MINUS)
sqtls = sqtls_raw.rename(columns={
    'Gene Symbol': 'gene',
    'Condition': 'condition',
    'P-value': 'pval',
    'Effect direction': 'effect_dir'
})
sqtls['pval_pro'] = sqtls['pval'].mask(sqtls['condition'] != 'Proliferative', 0)
sqtls['pval_qui'] = sqtls['pval'].mask(sqtls['condition'] != 'Quiescent', 0)
sqtls = sqtls.groupby(['gene_id', 'start', 'end']).agg(first).reset_index()
# sqtls = sqtls[sqtls['gene'].isin(['PROCR', 'DMAC2', 'NDUFAF6', 'NME7', 'PARP12'])]

# %%
def stem(accession):
    return accession.split('.')[0]
def similar(acc1, acc2):
    return stem(acc1) == stem(acc2)
gene_id_stems = {stem(gene_id) for gene_id in sqtls['gene_id']}
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
junc_colors = {
    (True, False): '#FEEFF4',
    (False, True): '#F0F0FD',
    (True, True): '#F8EEFD'
}

def sortkey(transcript: 'Transcript'):
    return -expression.average[transcript.accession], getattr(transcript, 'sqanti', SQANTI.OTHER)

def abundance(transcript: 'Transcript'):
    try:
        return expression.loc[transcript.accession]['average']
    except KeyError:
        return 0.0

def get_augmented_sqtl_record(row):
    chr, strand, gene_id, start, end, pval_pro, pval_qui, effect_dir = row
    if np.isnan(pval_pro):
        pval_pro = None
    if np.isnan(pval_qui):
        pval_qui = None
    try:
        gene = genes[stem(gene_id)]
    except KeyError:
        return None
    if strand is Strand.MINUS:
        start, end = end, start
    junc = Junction(start, end, chr, strand)
    junc_info = {
        'chr': chr,
        'strand': str(junc.strand),
        'donor': start,
        'acceptor': end,
        'gene': gene.name,
        'pval_pro': pval_pro,
        'pval_qui': pval_qui,
        'effect_dir': effect_dir
    }
    pb_transcripts = {tx for tx in gene.transcripts if isinstance(tx, PacBioTranscript) and all(orf.has_stop_codon for orf in tx.orfs)}
    if not pb_transcripts:
        return None
    lacking, containing = split_transcripts_on_junction_usage(junc, filter(attrgetter('orfs'), pb_transcripts))
    containing = sorted(containing, key=sortkey)
    lacking = sorted(lacking, key=sortkey)
    junc_info['containing'] = len(containing)
    junc_info['lacking'] = len(lacking)
    # if not containing or not lacking:
    #     return None
    
    gc_transcripts = sorted((tx for tx in gene.transcripts if isinstance(tx, GencodeTranscript)), key=attrgetter('appris'), reverse=True)
    anchor_tx = gc_transcripts[0]

    pairs = list(product(lacking, containing))
    # n_pairs = len(pairs)
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
    assert abs(sum(weights.values()) - 1) < 2**-8 or not weights

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

    freqs = {
        event: sum(int(event & pblock.flags == event)*weights[pblock.anchor.transcript, pblock.other.transcript] for pblock in pblocks)
        for event in AnnotationFlag.__members__.values() if event is not AnnotationFlag.NONE}
    junc_info.update(freqs)

    if not os.path.isdir(f'{output_dir}/{gene.name}'):
        os.mkdir(f'{output_dir}/{gene.name}')

    export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks)

    force_plotting = True
    fig_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
    if force_plotting or not os.path.isfile(fig_path):
        isoplot = IsoformPlot(
            gc_transcripts + [None] + containing + lacking + [None],
            columns = {
                'novelty': lambda tx: str(getattr(tx, 'sqanti', '')),
                'GENCODE': lambda tx: tx.gencode.name if getattr(tx, 'gencode', False) else '',
            }
        )
        with ExceptionLogger(f'Error plotting {gene.name} {junc}'):
            if not containing:  # adjust xlims to include unused junction if necessary
                a, b = junc.donor, junc.acceptor
                if b < a:
                    a, b = b, a
                SPACE = 30
                new_xlims = tuple(
                    interval for interval in ((x - SPACE, x + SPACE) for x in (a, b))
                    if not isoplot._subaxes.overlaps_range(*interval)
                )
                isoplot.xlims = isoplot.xlims + new_xlims
            gs = GridSpec(1, 3)
            isoplot.draw_all_isoforms(subplot_spec=gs[:-1])
            isoplot.draw_frameshifts(anchor=anchor_tx)
            isoplot.draw_region(type='line', color='k', linestyle='--', linewidth=1, track=len(gc_transcripts) + len(containing) + 0.5, start=gene.start, stop=gene.stop)
            junc_color = junc_colors[bool(pval_pro), bool(pval_qui)]
            isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor=junc_color)
            isoplot._handles['sQTL junction'] = mpatches.Patch(facecolor=junc_color)
            for pblock in pblocks:
                isoplot.draw_protein_block(pblock, alpha=weights[pblock.anchor.transcript, pblock.other.transcript]**0.5)
            isoplot.draw_features()
            # isoplot.draw_text(x=junc.acceptor+10, y=len(gc_transcripts), text=f'p-value: {pval}, {effect_dir}')
            isoplot.draw_legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
            heatmap_ax = isoplot.fig.add_subplot(gs[-1])
            all_accessions = [getattr(tx, 'accession', None) for tx in isoplot.transcripts]
            pb_accessions = [getattr(tx, 'accession', None) for tx in isoplot.transcripts if isinstance(tx, PacBioTranscript)]
            expression_sub = expression.loc[pb_accessions, :].iloc[:, list(range(0, 12, 2)) + list(range(1, 12, 2)) + [12, 13]]
            expression_sub = expression_sub.reindex(index=all_accessions, fill_value=0)
            expression_sub.iloc[-1, :] = expression_sub.sum(axis=0)
            heatmap = sns.heatmap(
                expression_sub,
                ax = heatmap_ax,
                xticklabels = True,
                yticklabels = False,
                linecolor = '#f0f0f0',
                linewidths = 1,
                cmap = sns.cubehelix_palette(as_cmap=True, light=0.995, dark=0.0),
                cbar_kws = {'label': 'CPM'},
                norm = LogNorm(vmin=0.1, vmax=10000.0),
                mask = expression_sub.applymap(lambda x: x == 0),
                annot = expression_sub,
                fmt = '0.0f',
            )
            heatmap.set_ylabel(None)
            isoplot._bax.set_title(f'p-value (proliferative): {pval_pro if pval_pro else ""}\np-value (quiescent): {pval_qui if pval_qui else ""}\n{effect_dir}', loc='left', size='x-small')
            isoplot.fig.set_size_inches(20, 0.8 + 0.4*len(gene.transcripts))
            plt.savefig(fig_path, facecolor='w', transparent=False, dpi=200, bbox_inches='tight')
            tqdm.write('\tsaved '+fig_path)
        plt.close(isoplot.fig)
    return junc_info

rows = [(row.chr, row.strand, row.gene_id, row.start, row.end, row.pval_pro, row.pval_qui, row.effect_dir) for row in sqtls.itertuples()]
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
