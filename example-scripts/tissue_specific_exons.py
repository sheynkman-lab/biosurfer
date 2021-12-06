#%%
from itertools import groupby
from operator import itemgetter

import matplotlib.pyplot as plt
import pandas as pd
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
from sqlalchemy import desc, func, select
from sqlalchemy.orm import contains_eager, joinedload, raiseload
from tqdm import tqdm

db = Database('gencode')
gtex = pd.read_csv('../data/gencode/b_gtex_isoform_medians_53tiss.tsv', sep='\t').set_index('TargetID').sort_values('Gene_Symbol')

#%%
try:
    frac_abundances = pd.read_csv('../output/ts-cassette-exons/frac_abundances.tsv', sep='\t')
except FileNotFoundError:
    tissues = gtex.columns[3:]
    total_abundances = gtex.drop(columns='Coord').groupby('Gene_Symbol', sort=False).agg('sum')

    def get_frac_abundance(col):
        gene_symbol = col[0]
        return col[3:].divide(total_abundances.loc[gene_symbol].replace(0.0, 1.0))

    frac_abundances = gtex.apply(get_frac_abundance, axis=1)
    frac_abundances.to_csv('../output/ts-cassette-exons/frac_abundances.tsv', sep='\t')

display(frac_abundances)

#%%
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
        where(Gene.chromosome_id == 'chr19').
        order_by(Transcript.gene_id)
    ).scalars())
    nrows = session.execute(
        select(func.count(Transcript.accession)).
        select_from(ORF).
        join(ORF.transcript).
        where(Transcript.gene_id.in_(gene_ids))
    ).scalars().first()

#%%
t = tqdm(None, desc='', total=nrows, unit='proteins', mininterval=0.2)
cassette_exon_pblocks = set()
intersects_domain = set()
intersects_idr = set()
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
                            cassette_exon_pblocks.add(pblock)
                            transcripts.append(other.transcript)
                            if pblock.category is PblockCategory.DELETION:
                                pblock_start = pblock.anchor_residues[0].position
                                pblock_stop = pblock.anchor_residues[-1].position
                            else:
                                pblock_start = pblock[0].anchor.position
                                pblock_stop = pblock_start
                            if any(pblock_start <= dom_stop and dom_start <= pblock_stop for dom_start, dom_stop in domains):
                                intersects_domain.add(pblock)
                            if any(pblock_start <= idr_stop and idr_start <= pblock_stop for idr_start, idr_stop in idrs):
                                intersects_idr.add(pblock)
            if len(transcripts) > 1:
                with ExceptionLogger(anchor):
                    isoplot = IsoformPlot(transcripts)
                    isoplot.draw_all_isoforms()
                    isoplot.draw_domains()
                    isoplot.draw_frameshifts()
                    isoplot.draw_legend()
                    fig_path = f'../output/ts-cassette-exons/{gene_name}.png'
                    plt.savefig(fig_path, facecolor='w', transparent=False, dpi=150, bbox_inches='tight')
                    tqdm.write('saved '+fig_path)
                    plt.close(isoplot.fig)
            t.update(len(proteins))
export_annotated_pblocks_to_tsv('../output/ts-cassette-exons/cassette_exon_pblocks.tsv', cassette_exon_pblocks)
export_annotated_pblocks_to_tsv('../output/ts-cassette-exons/intersects_domain.tsv', intersects_domain)
export_annotated_pblocks_to_tsv('../output/ts-cassette-exons/intersects_idr.tsv', intersects_idr)

# %%
