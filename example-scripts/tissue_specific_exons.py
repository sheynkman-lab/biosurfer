#%%
from itertools import chain, groupby
from operator import itemgetter

import matplotlib.pyplot as plt
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.alignments import Alignment, export_annotated_pblocks_to_tsv
from biosurfer.core.constants import AnnotationFlag, ProteinLevelAlignmentCategory as PblockCategory, FeatureType
from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import ORF, Gene, Protein, Transcript, Exon
from biosurfer.core.models.features import ProteinFeature
from biosurfer.plots.plotting import IsoformPlot
from more_itertools import chunked
from sqlalchemy import select, desc, func
from sqlalchemy.orm import contains_eager, joinedload, raiseload
from tqdm import tqdm

db = Database('gencode')

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
                    del_or_ins_pblocks = [pblock for pblock in aln.protein_blocks if pblock.category in {PblockCategory.DELETION, PblockCategory.INSERTION}]
                    if len(del_or_ins_pblocks) == 1:
                        pblock = del_or_ins_pblocks[0]
                        if pblock.event == 'SE' or pblock.event == 'IE':
                            cassette_exon_pblocks.add(pblock)
                            transcripts.append(other.transcript)
                            if pblock.category is PblockCategory.DELETION:
                                pblock_start = pblock.anchor_residues[0].position
                                pblock_stop = pblock.anchor_residues[-1].position
                            else:
                                pblock_start = pblock[0].anchor.position
                                pblock_stop = pblock_start
                            for dom_start, dom_stop in domains:
                                if pblock_start <= dom_stop and dom_start <= pblock_stop:
                                    intersects_domain.add(pblock)
                            for idr_start, idr_stop in idrs:
                                if pblock_start <= idr_stop and idr_start <= pblock_stop:
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
