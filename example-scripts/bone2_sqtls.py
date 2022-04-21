# %%
import csv
import multiprocessing as mp
import os
import re
import sys
from collections import Counter
from functools import reduce
from itertools import chain, groupby, product
from operator import attrgetter, itemgetter
from typing import TYPE_CHECKING, NamedTuple
from warnings import warn

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biosurfer.analysis.sqtl import (
    get_cblocks_attributed_to_transcript_event,
    get_transcript_events_associated_with_junction,
    split_transcripts_on_junction_usage)
from biosurfer.core.alignments import (CodonAlignment, ProteinAlignment,
                                       TranscriptAlignment)
from biosurfer.core.constants import OTHER_EXCLUSIVE, SQANTI, Strand
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (GencodeTranscript, Gene,
                                                Junction, PacBioTranscript,
                                                Transcript)
from biosurfer.plots.plotting import IsoformPlot, mpatches
from IPython.display import display
from matplotlib.colors import LogNorm
from matplotlib.gridspec import GridSpec
from more_itertools import first, only
from scipy.sparse import coo_matrix
from sqlalchemy.orm import sessionmaker
from sqlalchemy.sql.expression import and_, not_, or_, select
from tqdm import tqdm

if TYPE_CHECKING:
    from biosurfer.core.alignments import (CodonAlignmentBlock,
                                           ProteinAlignmentBlock)
    from biosurfer.core.models.features import ProteinFeature
    from biosurfer.core.splice_events import BasicTranscriptEvent

plt.switch_backend('agg')

data_dir = '../data/bone2'
output_dir = '../output/bone2'

db = Database('bone2')

# %%
print('Loading isoform expression table...')
timepoints = [f't{t}' for t in (0, 2, 4, 10)]
try:
    expression = pd.read_csv(f'{data_dir}/corrected_expression.tsv', sep='\t', index_col='rn')
except FileNotFoundError:
    expr_raw = pd.read_csv(f'{data_dir}/counts_threshold.tsv', sep='\t', index_col='rn', usecols=[0] + list(range(117-11, 117)))

    print('Loading transcript collapse mapping table...')
    collapsed = dict()
    with open(f'{data_dir}/hfobs_orf_refined.tsv') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            group = row['pb_accs'].split('|')
            base = row['base_acc']
            if base not in expr_raw.index:
                try:
                    base = first(acc for acc in group if acc in expr_raw.index)
                except ValueError:
                    continue
            for tx in group:
                collapsed[tx] = base

    print('Correcting isoform expression table...')
    expr_raw = expr_raw.filter(items=collapsed.keys(), axis=0)
    tx_to_idx = {tx: i for i, tx in enumerate(expr_raw.index)}
    # N = len(expr_raw.index)
    ijv = ((tx_to_idx[collapsed.get(tx, tx)], i, 1) for tx, i in tx_to_idx.items())
    row, col, data = zip(*ijv)
    transform = coo_matrix((data, (row, col))).tocsr()

    pattern = re.compile(r't\d+[abc]')
    with db.get_session() as session:
        all_pb_ids = set(acc for acc in session.execute(select(PacBioTranscript.accession)).scalars())
    expression = pd.DataFrame(data=transform.dot(expr_raw.to_numpy()), index=expr_raw.index, columns=expr_raw.columns)
    expression = (expression.
        filter(items=all_pb_ids, axis=0).
        rename(columns=lambda name: match[0] if (match := pattern.search(name)) else name)
    )
    for t in timepoints:
        expression[t] = expression[[col for col in expression.columns if col.startswith(t)]].mean(axis=1)
    expression['average'] = expression[timepoints].mean(axis=1)

    expression.to_csv(f'{data_dir}/corrected_expression.tsv', sep='\t')

# %%
print('Loading sQTL table...')
sqtls: pd.DataFrame = pd.read_csv(f'{data_dir}/sqtl_coloc.tsv', sep='\t', nrows=None)
sqtls[['chr', 'start', 'end', 'cluster', 'gene_id']] = sqtls['phenotype_id'].str.split(':', expand=True)
sqtls['gene_id_stem'] = sqtls['gene_id'].str.split('.').str.get(0)
sqtls[['start', 'end']] = sqtls[['start', 'end']].astype(int)
sqtls = sqtls.sample(frac=0.01, axis=0, random_state=329)
# sqtls = sqtls[sqtls['gene_id_stem'] == 'ENSG00000185340']

# %%
gene_id_mapper = dict()
with db.get_session() as session:
    for gene_id_stem in sqtls['gene_id_stem'].unique():
        gene_id = session.execute(
            select(Gene.accession).
            where(Gene.accession.startswith(gene_id_stem))
        ).scalar()
        if gene_id:
            gene_id_mapper[gene_id_stem] = gene_id

# db.project_feature_mappings(gene_ids=gene_id_mapper.values())

# %%
# def abundant_and_coding(transcript: 'Transcript'):
#     return transcript.accession in over3counts and len(transcript.orfs) > 0
junc_colors = {
    (True, False): '#FEEFF4',
    (False, True): '#F0F0FD',
    (True, True): '#F8EEFD'
}

def sortkey(transcript: 'Transcript'):
    return -abundance(transcript), getattr(transcript, 'sqanti', SQANTI.OTHER)

def abundance(transcript: 'Transcript') -> float:
    try:
        return expression.average[transcript.accession]
    except KeyError:
        return 0.0

class CBlockEventTuple(NamedTuple):
    cblock: 'CodonAlignmentBlock'
    event: 'BasicTranscriptEvent'
    anchor: 'Transcript'
    other: 'Transcript'

class BlockTuple(NamedTuple):
    pblock: 'ProteinAlignmentBlock'
    cblock: 'CodonAlignmentBlock'
    events: tuple['BasicTranscriptEvent', ...]
    anchor: 'Transcript'
    other: 'Transcript'

class FeatureTuple(NamedTuple):
    feature: 'ProteinFeature'
    cblock: 'CodonAlignmentBlock'
    affects_whole: bool
    anchor: 'Transcript'
    other: 'Transcript'

def get_augmented_sqtl_record(row):
    chr, gene_id, start, end, pval, slope, maf = row
    try:
        gene_id = gene_id_mapper[gene_id]
    except KeyError:
        return None
    
    with db.get_session() as session:
        gene = Gene.from_accession(session, gene_id)
        if not gene:
            return None
        start += 1
        end -= 1
        strand = gene.strand
        if strand is Strand.MINUS:
            start, end = end, start
        junc = Junction.from_coordinates(chr, strand, start, end)
        junc_info = {
            'chr': chr,
            'strand': str(junc.strand),
            'donor': start,
            'acceptor': end,
            'gene': gene.name,
            'pval': pval,
            'slope': slope,
            'maf': maf
        }
        pb_transcripts = {tx for tx in gene.transcripts if tx.accession in expression.index and any(orf.has_stop_codon for orf in tx.orfs)}
        if not pb_transcripts:
            return None
        lacking, containing = split_transcripts_on_junction_usage(junc, filter(attrgetter('orfs'), pb_transcripts))
        containing = sorted(containing, key=sortkey)
        lacking = sorted(lacking, key=sortkey)
        junc_info['containing'] = len(containing)
        junc_info['lacking'] = len(lacking)
        if not containing:
            return None
        
        gc_transcripts = sorted((tx for tx in gene.transcripts if isinstance(tx, GencodeTranscript)), key=attrgetter('appris'), reverse=True)
        anchor_tx = gc_transcripts[0]

        pairs = list(product(lacking, containing))
        
        abundance_denom = sum(abundance(tx) for tx in containing) * sum(abundance(tx) for tx in lacking)
        weights = {
            (tx1, tx2): abundance(tx1) * abundance(tx2) / abundance_denom
            for tx1, tx2 in pairs
        }
        if weights and abs(error := sum(weights.values()) - 1) > 2**-8:
            warn(f'Weights add up to 1{error:+.3e}')

        top_pair, junc_info['highest_pair_weight'] = max(weights.items(), key=itemgetter(1), default=(None, None))

        # calculate fraction of pairs where one isoform is NMD and other is not
        junc_info['NMD'] = None
        with ExceptionLogger(f'Error for {gene.name} {junc}'):
            junc_info['NMD'] = sum(
                float(tx1.primary_orf.nmd ^ tx2.primary_orf.nmd) * weights[tx1, tx2]
                for tx1, tx2 in pairs
            )

        tx_aln_to_events: dict['TranscriptAlignment', set['BasicTranscriptEvent']] = dict()
        for anchor, other in pairs:
            with ExceptionLogger(f'Error for {anchor} and {other}'):
                anchor.nucleotides, other.nucleotides
                tx_aln = TranscriptAlignment.from_transcripts(anchor, other)
                events = set(get_transcript_events_associated_with_junction(junc, tx_aln))
                tx_aln_to_events[tx_aln] = events
        cblock_event_tuples = [
            CBlockEventTuple(cblock, event, tx_aln.anchor, tx_aln.other)
            for tx_aln, events in tx_aln_to_events.items()
            for event in events
            for cblock in get_cblocks_attributed_to_transcript_event(event, CodonAlignment.from_proteins(tx_aln.anchor.protein, tx_aln.other.protein))
        ]
        cblock_event_tuples.sort(key=lambda cet: (cet.anchor.name, cet.other.name, cet.cblock, cet.event.start, cet.event.stop))
        block_tuples = [
            BlockTuple(
                ProteinAlignment.from_proteins(anchor.protein, other.protein).cblock_to_pblock[cblock],
                cblock,
                tuple(cet.event for cet in group),
                anchor, other
            )
            for (anchor, other, cblock), group in groupby(cblock_event_tuples, key=lambda cet: (cet.anchor, cet.other, cet.cblock))
        ]

        feature_tuples = list(chain(
            (
                FeatureTuple(
                    feature,
                    bt.cblock,
                    bt.cblock.anchor_range.start <= feature.protein_start - 1 and feature.protein_stop <= bt.cblock.anchor_range.stop,
                    bt.anchor,
                    bt.other
                )
                for bt in block_tuples if (bt.anchor, bt.other) == top_pair
                for feature in bt.anchor.protein.features if feature.protein_start - 1 < bt.cblock.anchor_range.stop and bt.cblock.anchor_range.start < feature.protein_stop
            ),
            (
                FeatureTuple(
                    feature,
                    bt.cblock,
                    bt.cblock.other_range.start <= feature.protein_start - 1 and feature.protein_stop <= bt.cblock.other_range.stop,
                    bt.anchor,
                    bt.other
                )
                for bt in block_tuples if (bt.anchor, bt.other) == top_pair
                for feature in bt.other.protein.features if feature.protein_start - 1 < bt.cblock.other_range.stop and bt.cblock.other_range.start < feature.protein_stop
            )
        ))

        delta_aa_weights = tuple(zip(*(
            (sum(bt.cblock.delta_length for bt in group), weights[pair])
            for pair, group in groupby(block_tuples, key=lambda bt: (bt.anchor, bt.other))
        )))
        if delta_aa_weights:
            junc_info['avg_delta_aa'] = np.average(delta_aa_weights[0], weights=delta_aa_weights[1])
        else:
            junc_info['avg_delta_aa'] = None

        # calculate event frequencies
        event_freqs = reduce(
            Counter.__add__,
            (
                Counter({k: v*weight for k, v in counts.items()})
                for counts, weight in (
                    (Counter(events), weights[tx_aln.anchor, tx_aln.other])
                    for tx_aln, events in tx_aln_to_events.items()
                )
            ),
            Counter()
        )
        junc_info['most_common_event'] = ', '.join(
            f'{type(event).__name__}(start={event.start.coordinate}, stop={event.stop.coordinate})'
            for event, freq in event_freqs.items() if freq == max(event_freqs.values())
        )
        junc_info['highest_event_freq'] = max(event_freqs.values()) if event_freqs else None

        junc_info['affected_features'] = {ft.feature.feature.name for ft in feature_tuples}

        if not os.path.isdir(f'{output_dir}/{gene.name}'):
            os.mkdir(f'{output_dir}/{gene.name}')

        blocks_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor.coordinate}_{junc.acceptor.coordinate}_blocks.tsv'
        with open(blocks_path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['protein_alignment_block', 'codon_alignment_block', 'delta_aa', 'transcript_events', 'anchor', 'other', 'anchor_cpm', 'other_cpm', 'pair_weight'])
            writer.writerows(
                (
                    bt.pblock,
                    bt.cblock,
                    bt.cblock.delta_length,
                    ', '.join(f'{type(event).__name__}(start={event.start.coordinate}, stop={event.stop.coordinate})' for event in bt.events),
                    bt.anchor, bt.other,
                    abundance(bt.anchor), abundance(bt.other),
                    weights[bt.anchor, bt.other]
                )
                for bt in block_tuples
            )
        
        features_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor.coordinate}_{junc.acceptor.coordinate}_features.tsv'
        with open(features_path, 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerow(['affected_feature', 'codon_alignment_block', 'affects_whole_feature', 'anchor', 'other'])
            writer.writerows(feature_tuples)

        force_plotting = False
        fig_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor.coordinate}_{junc.acceptor.coordinate}.png'

        if force_plotting or not os.path.isfile(fig_path):
            isoplot = IsoformPlot(
                gc_transcripts + [None] + containing + [None] + lacking + [None, None],
                columns = {
                    'novelty': lambda tx: str(getattr(tx, 'sqanti', '')),
                    'GENCODE': lambda tx: tx.gencode.name if getattr(tx, 'gencode', False) else '',
                }
            )
            with ExceptionLogger(f'Error plotting {gene.name} {junc}'):
                # if not containing:  # adjust xlims to include unused junction if necessary
                #     c, b = junc.donor, junc.acceptor
                #     if b < c:
                #         c, b = b, c
                #     SPACE = 30
                #     new_xlims = tuple(
                #         interval for interval in ((x - SPACE, x + SPACE) for x in (c, b))
                #         if not isoplot._subaxes.overlaps_range(*interval)
                #     )
                #     isoplot.xlims = isoplot.xlims + new_xlims
                gs = GridSpec(1, 5)
                isoplot.draw_all_isoforms(subplot_spec=gs[:-1])
                isoplot.draw_frameshifts(anchor=anchor_tx)
                isoplot.draw_region(type='line', color='k', linewidth=1, track=len(gc_transcripts), start=gene.start, stop=gene.stop)
                isoplot.draw_region(type='line', color='k', linestyle='--', linewidth=1, track=len(gc_transcripts) + len(containing) + 1, start=gene.start, stop=gene.stop)
                junc_color = '#eeeeee'
                isoplot.draw_background_rect(start=junc.donor.coordinate, stop=junc.acceptor.coordinate, facecolor=junc_color)
                isoplot._handles['sQTL junction'] = mpatches.Patch(facecolor=junc_color)
                for (anchor, other), group in groupby(block_tuples, key=lambda bt: (bt.anchor, bt.other)):
                    isoplot.draw_protein_alignment_blocks((pt[0] for pt in group), anchor.protein, other.protein, alpha = weights[anchor, other]*0.75)
                isoplot.draw_features()
                isoplot.draw_legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
                heatmap_ax = isoplot.fig.add_subplot(gs[-1])
                all_accessions = [getattr(tx, 'accession', None) for tx in isoplot.transcripts]
                pb_accessions = [getattr(tx, 'accession', None) for tx in isoplot.transcripts if isinstance(tx, PacBioTranscript)]
                expression_sub = expression.loc[pb_accessions, timepoints + ['average']]
                expression_sub = expression_sub.reindex(index=all_accessions, fill_value=0)
                expression_sub.iloc[-1, :] = expression_sub.sum(axis=0)
                # c = len(gc_transcripts) + len(containing) + 2
                # expression_sub.iloc[c, :] = expression_sub.iloc[:c, :].sum(axis=0)
                # expression_sub.iloc[-3, :] = expression_sub.iloc[c+1:-3, :].sum(axis=0)
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
                isoplot._bax.set_title(f'p-value: {pval:.3e}\nslope: {slope:.3f}\nMAF: {maf:.3f}', loc='left', size='x-small')
                isoplot.fig.set_size_inches(20, 0.8 + 0.4*len(isoplot.transcripts))
                plt.savefig(fig_path, facecolor='w', transparent=False, dpi=200, bbox_inches='tight')
                # tqdm.write('\tsaved '+fig_path)
            plt.close(isoplot.fig)
    return junc_info

rows = [(row.chr, row.gene_id_stem, row.start, row.end, row.pval_nominal, row.slope, row.maf) for row in sqtls.itertuples()]
records = []
with mp.Pool() as p:
    with (
    tqdm(desc='Analyzing sQTL junctions', total=sqtls.shape[0], file=sys.stdout, unit='junctions') as t,
    tqdm(desc='Annotated junctions', file=sys.stdout, unit='junctions', miniters=1) as t2):
        for result in p.imap_unordered(get_augmented_sqtl_record, rows):
            t.update()
            if result:
                records.append(result)
                t2.update()

sqtls_augmented = pd.DataFrame.from_records(records)
# display(sqtls_augmented)
# print(f'Annotated {sqtls_augmented.shape[0]} sQTL junctions')
sqtls_augmented.to_csv(f'{output_dir}/coloc_sqtls_annotated.tsv', sep='\t', index=False)

# %%
