#%%
from itertools import chain, groupby
from operator import itemgetter

from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.alignments import Alignment, export_annotated_pblocks_to_tsv
from biosurfer.core.constants import AnnotationFlag, ProteinLevelAlignmentCategory as PblockCategory
from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import ORF, Gene, Protein, Transcript, Exon
from biosurfer.core.models.features import ProteinFeature
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
cassette_exon_pblocks = []
with db.get_session() as session:
    for gene_chunk in chunked(gene_ids, 200):
        rows = session.execute(
            q.where(Transcript.gene_id.in_(gene_chunk))
        ).unique()
        rows_by_gene = groupby(rows, key=itemgetter(1))
        for gene, group in rows_by_gene:
            proteins = [row[0] for row in group]
            anchor = proteins[0]
            for other in proteins[1:]:
                with ExceptionLogger(f'{anchor}|{other}'):
                    aln = Alignment(anchor, other)
                    del_or_ins_pblocks = [pblock for pblock in aln.protein_blocks if pblock.category in {PblockCategory.DELETION, PblockCategory.INSERTION}]
                    if len(del_or_ins_pblocks) == 1:
                        pblock = del_or_ins_pblocks[0]
                        if pblock.event == 'SE' or pblock.event == 'IE':
                            cassette_exon_pblocks.append(pblock)
            t.update(len(proteins))
export_annotated_pblocks_to_tsv('../output/cassette_exon_pblocks.tsv', cassette_exon_pblocks)

# %%
