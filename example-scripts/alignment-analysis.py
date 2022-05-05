import multiprocessing as mp
from itertools import chain, groupby
from operator import attrgetter
from typing import Iterable

import pandas as pd
from biosurfer.core.alignments import ProteinAlignment, SeqAlignCat
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import Gene, Transcript, ORF
from more_itertools import only, first, last
from sqlalchemy import select
from sqlalchemy.orm import contains_eager, joinedload

db = Database('gencode')

def process_gene(transcripts: Iterable['Transcript']):
    out = []
    transcripts = sorted(transcripts, key=attrgetter('appris'), reverse=True)
    principal = transcripts[0]
    if principal.appris is not APPRIS.PRINCIPAL:
        return out    
    for alternative in transcripts[1:]:
        pr_aln = ProteinAlignment.from_proteins(principal.protein, alternative.protein)
        pblocks = pr_aln.blocks
        internal = [pblock for pblock in pblocks if pblock.category is not SeqAlignCat.MATCH]
        nterm = only(pblock for pblock in internal if pblock is first(pblocks, None))
        cterm = only(pblock for pblock in internal if pblock is last(pblocks, None))
        internal = [pblock for pblock in internal if not (pblock is nterm or pblock is cterm)]
        out.append({
            'anchor': principal.name,
            'other': alternative.name,
            'num_internal': len(internal),
            'has_nterm': bool(nterm),
            'has_cterm': bool(cterm)
        })
    return out

genes = ['PAQR6']
with db.get_session() as session:
    transcripts = session.execute(
        select(Transcript).
        select_from(ORF).
        join(ORF.transcript).
        join(Transcript.gene).
        where(Gene.name.in_(genes))
    ).scalars().all()
    records = list(chain.from_iterable(map(process_gene, (txs for _, txs in groupby(transcripts, key=attrgetter('gene_id'))))))
df = pd.DataFrame.from_records(records)
