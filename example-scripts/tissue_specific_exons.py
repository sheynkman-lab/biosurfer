#%%
import csv
import os
from itertools import groupby
from operator import itemgetter

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biosurfer.core.alignments import (Alignment,
                                       export_annotated_pblocks_to_tsv)
from biosurfer.core.constants import FeatureType
from biosurfer.core.constants import \
    ProteinLevelAlignmentCategory as PblockCategory
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (ORF, Exon, Gene, Protein,
                                                Transcript)
from biosurfer.core.models.features import ProteinFeature
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from more_itertools import chunked
from scipy.stats import chi2_contingency
from sklearn.linear_model import LogisticRegression
from sqlalchemy import desc, func, select
from sqlalchemy.orm import contains_eager, joinedload, raiseload
from tqdm import tqdm

db = Database('gencode')
gtex = pd.read_csv('../data/gencode/b_gtex_isoform_medians_53tiss.tsv', sep='\t').set_index('TargetID').sort_values('Gene_Symbol')

def strip_version(accession):
    return '.'.join(accession.split('.')[:-1])

gtex.index = gtex.index.map(strip_version)

#%%
try:
    frac_abundances = pd.read_csv('../output/ts-cassette-exons/frac_abundances.tsv', sep='\t', index_col='TargetID')
except FileNotFoundError:
    tissues = gtex.columns[3:]
    total_abundances = gtex.drop(columns='Coord').groupby('Gene_Symbol', sort=False).agg('sum')

    def get_frac_abundance(col):
        gene_symbol = col[0]
        return col[3:].divide(total_abundances.loc[gene_symbol].replace(0.0, 1.0))

    frac_abundances = gtex.apply(get_frac_abundance, axis=1)
    frac_abundances.to_csv('../output/ts-cassette-exons/frac_abundances.tsv', sep='\t')

frac_abundances.index = frac_abundances.index.map(strip_version)
enst_to_ensg = dict(gtex['Gene_Symbol'])

display(frac_abundances)

#%%
genes = frac_abundances.index.map(enst_to_ensg.get)
normalized_frac_abundances = frac_abundances.divide(frac_abundances.max(axis=1).replace(0.0, 1.0), axis=0)
display(normalized_frac_abundances)

tau = (1 - normalized_frac_abundances).sum(axis=1) / (len(normalized_frac_abundances.columns) - 1)
tau[tau > 1.0] = np.nan
display(tau)

# nbins = 4
# enst_to_bin, bins = pd.cut(tau, nbins, retbins=True)
# enst_to_bin.dropna(inplace=True)
# bin_counts = enst_to_bin.value_counts()

# contingency_domains = pd.DataFrame(
#     {'non-intersecting': [0 for _ in range(nbins)], 'intersecting': [0 for _ in range(nbins)]},
#     index = bin_counts.index
# )
# contingency_domains.index.name = 'tau'
# contingency_idrs = contingency_domains.copy()

#%%
temp_file = '../output/ts-cassette-exons/cassette_exon_ensts.tsv'
if not os.path.isfile(temp_file):
    with db.get_session() as session:
        protein_tx = contains_eager(Protein.orf).contains_eager(ORF.transcript)
        q = (
            select(Protein, Transcript.gene_id).
            join(Protein.orf).
            join(ORF.transcript).
            order_by(Transcript.gene_id).
            order_by(desc(Transcript.__table__.c.appris)).
            order_by(desc(ORF.length)).
            order_by(Transcript.name).
            options(
                protein_tx.joinedload(Transcript.gene),
                protein_tx.joinedload(Transcript.orfs),
                protein_tx.joinedload(Transcript.exons).joinedload(Exon.transcript),
                contains_eager(Protein.orf).joinedload(ORF.protein),
                joinedload(Protein.features).joinedload(ProteinFeature.protein),
                joinedload(Protein.features).joinedload(ProteinFeature.feature),
                raiseload('*')
            )
        )
        tqdm.write(str(q))
        gene_ids = list(session.execute(
            select(Transcript.gene_id).distinct().
            select_from(ORF).
            join(ORF.transcript).
            join(Transcript.gene).
            # where(Gene.chromosome_id == 'chr19').
            order_by(Transcript.gene_id)
        ).scalars())
        nrows = session.execute(
            select(func.count(Transcript.accession)).
            select_from(ORF).
            join(ORF.transcript).
            where(Transcript.gene_id.in_(gene_ids))
        ).scalars().first()
        
    FIELDNAMES = ('TargetID', 'affects_domain', 'affects_idr')
    with open(temp_file, 'w') as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t')
        writer.writeheader()

    FLUSH_SIZE = 1000
    records = []
    t = tqdm(None, desc='', total=nrows, unit='proteins', mininterval=0.2)
    with db.get_session() as session:
        for gene_chunk in chunked(gene_ids, 200):
            rows = session.execute(
                q.where(Transcript.gene_id.in_(gene_chunk))
            ).unique()
            rows_by_gene = groupby(rows, key=itemgetter(1))
            for gene_id, group in rows_by_gene:
                proteins = [row[0] for row in group]
                anchor = proteins[0]
                gene_name = anchor.transcript.gene.name
                domains = tuple((f.protein_start, f.protein_stop) for f in anchor.features if f.type is FeatureType.DOMAIN)
                idrs = tuple((f.protein_start, f.protein_stop) for f in anchor.features if f.type is FeatureType.IDR)
                transcripts = [anchor.transcript]
                for other in proteins[1:]:
                    with ExceptionLogger(f'{anchor}|{other}'):
                        aln = Alignment(anchor, other)
                        nonmatch_pblocks = [pblock for pblock in aln.protein_blocks if pblock.category is not PblockCategory.MATCH]
                        if len(nonmatch_pblocks) == 1 and (pblock := nonmatch_pblocks[0]).category is not PblockCategory.SUBSTITUTION:
                            if pblock.event == 'SE' and ' to ' not in pblock.annotation or pblock.event == 'IE' and 'exons' not in pblock.annotation:
                                stripped_accession = strip_version(other.transcript.accession)
                                record = {'TargetID': stripped_accession}
                                transcripts.append(other.transcript)
                                if pblock.category is PblockCategory.DELETION:
                                    pblock_start = pblock.anchor_residues[0].position
                                    pblock_stop = pblock.anchor_residues[-1].position
                                else:
                                    pblock_start = pblock[0].anchor.position
                                    pblock_stop = pblock_start
                                record['affects_domain'] = any(pblock_start <= dom_stop and dom_start <= pblock_stop for dom_start, dom_stop in domains)
                                record['affects_idr'] = any(pblock_start <= idr_stop and idr_start <= pblock_stop for idr_start, idr_stop in idrs)
                                records.append(record)
                if len(records) >= FLUSH_SIZE:
                    with open(temp_file, 'a') as f:
                        writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t')
                        writer.writerows(records)
                        records[:] = []
    #             fig_path = f'../output/ts-cassette-exons/{gene_name}.png'
    #             if len(transcripts) > 1 and not os.path.isfile(fig_path):
    #                 with ExceptionLogger(anchor):
    #                     isoplot = IsoformPlot(transcripts)
    #                     isoplot.draw_all_isoforms()
    #                     isoplot.draw_domains()
    #                     isoplot.draw_frameshifts()
    #                     isoplot.draw_legend()
    #                     plt.savefig(fig_path, facecolor='w', transparent=False, dpi=150, bbox_inches='tight')
    #                     tqdm.write('saved '+fig_path)
    #                     plt.close(isoplot.fig)
                t.update(len(proteins))
            
        with open(temp_file, 'a') as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES, delimiter='\t')
            writer.writerows(records)

cassette_exon_txs = pd.read_csv(temp_file, sep='\t').set_index('TargetID')

# %%
cassette_exon_txs['tau'] = tau
cassette_exon_txs['category'] = pd.cut(cassette_exon_txs['tau'], bins=tau.quantile([0.0, 0.2, 0.8, 1.0]), labels=['constitutive', 'other', 'tissue-specific'])
cassette_exon_txs.dropna(inplace=True)
display(cassette_exon_txs)

# %%
contingency_domain = pd.crosstab(cassette_exon_txs['category'], cassette_exon_txs['affects_domain'])
display(contingency_domain)
display(chi2_contingency(contingency_domain))

contingency_idr = pd.crosstab(cassette_exon_txs['category'], cassette_exon_txs['affects_idr'])
display(contingency_idr)
display(chi2_contingency(contingency_idr))

# %%
plt.figure(1)
sns.violinplot(x='affects_domain', y='tau', data=cassette_exon_txs, cut=0)
plt.savefig('../output/ts-cassette-exons/domain_boxplot.png')

plt.figure(2)
sns.violinplot(x='affects_idr', y='tau', data=cassette_exon_txs, cut=0)
plt.savefig('../output/ts-cassette-exons/idr_boxplot.png')

# %%
