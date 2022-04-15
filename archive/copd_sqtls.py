# %%
import csv
import os
import traceback
from itertools import chain, groupby, product
from operator import attrgetter, itemgetter
from statistics import median
from warnings import filterwarnings

import matplotlib.pyplot as plt
import pandas as pd
from biosurfer.analysis.sqtl import (get_event_counts,
                                     get_pblocks_attributed_to_junction,
                                     junction_has_drastic_effect_in_pair,
                                     split_transcripts_on_junction_usage)
from biosurfer.core.alignments import (Alignment,
                                       export_annotated_pblocks_to_tsv)
from biosurfer.core.constants import Strand
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger, IntervalTree
from biosurfer.core.models.biomolecules import (ORF, Exon, Gene, Junction,
                                                Protein, Transcript)
from biosurfer.core.models.features import Feature, ProteinFeature
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning
from more_itertools import chunked
from sqlalchemy import select
from sqlalchemy.orm import joinedload, raiseload
from sqlalchemy.sql.expression import and_, func, or_
from tqdm import tqdm

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

db = Database('gencode')
session = db.get_session()

data_dir = '../data/copd'
output_dir = '../output/copd'
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

# %%
tqdm.write('Reading phenotype IDs...')
raw = pd.read_csv(f'{data_dir}/LTRC_sqtl_GWAS_0.1FDR_result.tsv', sep='\t', usecols=['phenotype_id', 'pval_nominal'], nrows=None)
sqtl_pvals = raw.groupby('phenotype_id').agg({'pval_nominal': 'min'}).nsmallest(n=500, columns='pval_nominal')
sqtls_raw = list(sqtl_pvals.index.unique())
sqtls = dict()
for id in sqtls_raw:
    chr, start, stop, cluster = id.split(':')
    sqtls[id] = {
        'chr': 'chr'+chr,
        'start': int(start),
        'stop': int(stop),
        'cluster': cluster.split('_')[1],
        'strand': cluster[-1]
    }

tqdm.write('Loading gene ranges...')
subq = (
    select(Gene.accession, Gene.chromosome_id, Exon.start, Exon.stop).
    select_from(ORF).
    join(ORF.transcript).
    join(Transcript.gene).
    join(Transcript.exons).
    where(Gene.chromosome_id.in_(sqtl['chr'] for sqtl in sqtls.values())).
    subquery()
)
gene_coords = session.execute(
    select(subq.c.accession, subq.c.chromosome_id, func.min(subq.c.start), func.max(subq.c.stop)).
    group_by(subq.c.accession).
    order_by(subq.c.chromosome_id)
)
chr_coord_to_gene = dict()
for chr, group in groupby(gene_coords, key=itemgetter(1)):
    rows = list(group)
    chr_coord_to_gene[chr] = IntervalTree.from_tuples((int(row[2]), int(row[3]), row[0]) for row in rows)

sqtl_to_genes = {id: [i[-1] for i in chr_coord_to_gene[sqtl['chr']][sqtl['start']:sqtl['stop']]] for id, sqtl in sqtls.items()}
gene_to_sqtls = dict()
for sqtl, genes in sqtl_to_genes.items():
    for gene in genes:
        gene_to_sqtls.setdefault(gene, []).append(sqtl)

# db.project_feature_mappings(gene_to_sqtls.keys())

# %%
records = []
t = tqdm(None, total=sum(len(gene_sqtls) for gene_sqtls in gene_to_sqtls.values()))
gene_tx = joinedload(Gene.transcripts)
gene_protein = gene_tx.joinedload(Transcript.orfs).joinedload(ORF.protein)
for gene_chunk in chunked(gene_to_sqtls.keys(), 100):
    t.desc = 'Loading objects from database'
    stmt = (
        select(Gene).
        where(Gene.accession.in_(gene_chunk)).
        options(
            gene_tx.joinedload(Transcript.gene),
            gene_tx.joinedload(Transcript.exons).joinedload(Exon.transcript),
            gene_tx.joinedload(Transcript.orfs).joinedload(ORF.transcript),
            gene_protein.joinedload(Protein.orf),
            gene_protein.joinedload(Protein.features).joinedload(ProteinFeature.protein),
            gene_protein.joinedload(Protein.features).joinedload(ProteinFeature.feature),
            raiseload('*')
        )
    )
    gene_dict = {gene.accession: gene for gene in session.execute(stmt).unique().scalars()}
    t.desc = 'Analyzing sQTLs'
    for gene_id, gene in gene_dict.items():
        if not os.path.isdir(f'{output_dir}/{gene.name}'):
            os.mkdir(f'{output_dir}/{gene.name}')
        for id in gene_to_sqtls[gene_id]:
            t.update()
            sqtl = sqtls[id]
            chr, donor, acceptor = sqtl['chr'], sqtl['start'], sqtl['stop']
            strand = Strand.from_symbol(sqtl['strand'])
            if strand is Strand.MINUS:
                donor, acceptor = acceptor, donor
            junc = Junction(donor, acceptor, chr, strand)
            coding_transcripts = sorted(
                (transcript for transcript in gene.transcripts if transcript.sequence and transcript.orfs),
                key = attrgetter('appris'),
                reverse = True
            )

            junc_info = {
                'chr': str(junc.chromosome),
                'strand': str(junc.strand),
                'donor': junc.donor,
                'acceptor': junc.acceptor,
                'gene': gene.name,
                'pairs': 0,
                'min_pval': float(sqtl_pvals['pval_nominal'][id])
            }

            not_using, using = split_transcripts_on_junction_usage(junc, coding_transcripts)
            using = list(using)
            not_using = list(not_using)

            if not all((using, not_using)):
                continue
            
            pairs = list(product(not_using, using))
            n_pairs = len(pairs)
            alns = [Alignment(tx1.protein, tx2.protein) for tx1, tx2 in pairs]
            junc_info['pairs'] = n_pairs

            if n_pairs > 0:
                pblocks = get_pblocks_attributed_to_junction(junc, alns)

                # calculate fraction of pairs where junction-lacking isoform is not NMD but junction-using isoform is NMD
                junc_info['NMD'] = sum(tx2.primary_orf.nmd and not tx1.primary_orf.nmd for tx1, tx2 in pairs) / n_pairs

                # calculate median change in sequence length for junction-related pblocks
                junc_info['median_delta_length'] = median(pblock.delta_length for pblock in pblocks) if pblocks else 0.0

                # classify knockdown vs. alteration for each pair
                # threshold = -not_using[0].protein.length * 2 // 5
                pair_pblocks = groupby(pblocks, key=attrgetter('anchor', 'other'))
                knockdown_pair_count = sum(junction_has_drastic_effect_in_pair(pblocks=list(pblock_group)) for _, pblock_group in pair_pblocks)
                junc_info['knockdown_frequency'] = knockdown_pair_count / n_pairs

                counts = get_event_counts(pblocks)
                freqs = {event.name: count/n_pairs for event, count in counts.items()}
                junc_info.update(freqs)
            else:
                pblocks = []

            records.append(junc_info)

            export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks)

            force_plotting = True
            fig_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
            if force_plotting or not os.path.isfile(fig_path):
                isoplot = IsoformPlot(using + not_using, columns={'APPRIS': attrgetter('appris.name')})
                try:
                    isoplot.draw_all_isoforms()
                    isoplot.draw_frameshifts()
                    isoplot.draw_features()
                    isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#f0f0f0')
                    isoplot.draw_region(type='line', color='k', linestyle='--', linewidth=1, track=len(using) - 0.5, start=gene.start, stop=gene.stop)
                    for pblock in pblocks:
                        isoplot.draw_protein_block(pblock, alpha=n_pairs**-0.5)
                    isoplot.draw_legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
                    isoplot.fig.set_size_inches(9, 0.3*len(gene.transcripts))
                    plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
                    tqdm.write('\tsaved '+fig_path)
                except Exception as e:
                    tqdm.write(f'\tcould not plot {gene.name}_{junc.donor}_{junc.acceptor}:')
                    tqdm.write(traceback.format_exc())
                plt.close(isoplot.fig)
        if not os.listdir(f'{output_dir}/{gene.name}'):
            os.rmdir(f'{output_dir}/{gene.name}')

sqtls_augmented = pd.DataFrame.from_records(records)
# display(sqtls_augmented)
sqtls_augmented.to_csv(f'{output_dir}/copd_sqtls_annotated.tsv', sep='\t', index=False)
