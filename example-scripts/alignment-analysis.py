# %%
import multiprocessing as mp
from itertools import chain, groupby
from operator import attrgetter
import sys
from typing import Iterable

import pandas as pd
from biosurfer.core.alignments import (CodonAlignment, ProteinAlignment,
                                       SeqAlignCat, TranscriptAlignment)
from biosurfer.core.constants import APPRIS
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (ORF, GencodeTranscript, Gene,
                                                Protein, Transcript)
from biosurfer.core.splice_events import get_event_code
from IPython.display import display
from sqlalchemy import select
from sqlalchemy.orm import contains_eager, joinedload
from tqdm import tqdm

db = Database('gencode')

# %%
gene_to_transcripts: dict[str, list[str]] = dict()
with db.get_session() as session:
    rows = session.execute(
            select(Gene.name, Transcript.name).
            select_from(Protein).
            join(Protein.orf).
            join(ORF.transcript).
            join(Transcript.gene).
            where(Gene.chromosome_id == 'chr19')
    ).all()
    for gene_name, tx_name in rows:
        gene_to_transcripts.setdefault(gene_name, []).append(tx_name)

def process_gene(gene_name: str):
    out = []
    with db.get_session() as session:
        transcripts = list(GencodeTranscript.from_names(session, gene_to_transcripts[gene_name]).values())
        transcripts.sort(key=attrgetter('appris'), reverse=True)
        principal = transcripts[0]
        if principal.appris is APPRIS.PRINCIPAL:
            for alternative in transcripts[1:]:
                pblocks = ()
                with ExceptionLogger(f'{principal}, {alternative}'):
                    tx_aln = TranscriptAlignment.from_transcripts(principal, alternative)
                    cd_aln = CodonAlignment.from_proteins(principal.protein, alternative.protein)
                    pr_aln = ProteinAlignment.from_proteins(principal.protein, alternative.protein)
                    pblocks = pr_aln.blocks                    
                for pblock in pblocks:
                    if pblock.category is SeqAlignCat.MATCH:
                        continue
                    for cblock in pr_aln.pblock_to_cblocks[pblock]:
                        tblock = cd_aln.cblock_to_tblock.get(cblock)
                        events = tx_aln.block_to_events.get(tblock, ())
                        out.append({
                            'anchor': principal.name,
                            'other': alternative.name,
                            'nterm': pblock is pblocks[0],
                            'cterm': pblock is pblocks[-1],
                            'pblock': pblock,
                            'cblock': cblock,
                            'tblock': tblock,
                            'events': get_event_code(events)
                        })
    return out

records = []
with mp.Pool() as p:
    t = tqdm(desc='Processing genes', total=len(gene_to_transcripts), unit='gene', file=sys.stdout)
    for result in p.imap_unordered(process_gene, gene_to_transcripts.keys()):
        records.extend(result)
        t.update()
df = pd.DataFrame.from_records(records)
df.to_csv('../output/alignment-analysis.tsv', sep='\t', index=False)

# %%
