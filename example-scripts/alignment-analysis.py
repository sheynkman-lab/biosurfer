# %%
import multiprocessing as mp
import os
import sys
from itertools import chain, groupby, starmap
from operator import attrgetter

import matplotlib.pyplot as plt
import pandas as pd
from biosurfer.core.alignments import (CodonAlignment, ProteinAlignment,
                                       SeqAlignCat, TranscriptAlignment)
from biosurfer.core.constants import APPRIS, CTerminalChange, NTerminalChange
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (ORF, GencodeTranscript, Gene,
                                                Protein, Transcript)
from biosurfer.core.splice_events import SpliceEvent, get_event_code
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from more_itertools import first, one
from sqlalchemy import select
from tqdm import tqdm

db = Database('gencode')

# %%
print('Loading gene and transcript names...')
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

df_path = '../output/alignment-analysis.tsv'
try:
    df = pd.read_csv(df_path, sep='\t')
except:
    def process_gene(gene_name: str):
        out = []
        with db.get_session() as session:
            transcripts: list['GencodeTranscript'] = list(GencodeTranscript.from_names(session, gene_to_transcripts[gene_name]).values())
            transcripts.sort(key=attrgetter('appris'), reverse=True)
            principal = transcripts[0]
            if principal.appris is APPRIS.PRINCIPAL:
                anchor_start_codon = principal.get_genome_coord_from_transcript_coord(principal.primary_orf.transcript_start - 1)
                anchor_stop_codon = principal.get_genome_coord_from_transcript_coord(principal.primary_orf.transcript_stop - 1)
                for alternative in transcripts[1:]:
                    other_start_codon = alternative.get_genome_coord_from_transcript_coord(alternative.primary_orf.transcript_start - 1)
                    other_stop_codon = alternative.get_genome_coord_from_transcript_coord(alternative.primary_orf.transcript_stop - 1)
                    anchor_starts_upstream = anchor_start_codon <= other_start_codon
                    anchor_stops_upstream = anchor_stop_codon <= other_stop_codon
                    pblocks = ()
                    with ExceptionLogger(info=f'{principal}, {alternative}', callback=lambda x, y, z: transcripts.remove(alternative)):
                        tx_aln = TranscriptAlignment.from_transcripts(principal, alternative)
                        cd_aln = CodonAlignment.from_proteins(principal.protein, alternative.protein)
                        pr_aln = ProteinAlignment.from_proteins(principal.protein, alternative.protein)
                        pblocks = pr_aln.blocks
                        anchor_start_cblock = one(cd_aln.anchor_blocks.at(0)).data
                        other_start_cblock = one(cd_aln.other_blocks.at(0)).data
                        anchor_stop_cblock = one(cd_aln.anchor_blocks.at(principal.protein.length - 1)).data
                        other_stop_cblock = one(cd_aln.other_blocks.at(alternative.protein.length - 1)).data
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
                                'pblock': str(pblock),
                                'cblock': str(cblock),
                                'tblock': str(tblock),
                                'events': get_event_code(events),
                                'complex_splicing': any(set(events).intersection(compound_event.members) for compound_event in tx_aln.events if isinstance(compound_event, SpliceEvent) and len(compound_event.members) > 1),
                                'affects_upstream_start': anchor_starts_upstream and cblock is anchor_start_cblock or not anchor_starts_upstream and cblock is other_start_cblock,
                                'affects_downstream_start': anchor_starts_upstream and cblock is other_start_cblock or not anchor_starts_upstream and cblock is anchor_start_cblock,
                                'affects_upstream_stop': anchor_stops_upstream and cblock is anchor_stop_cblock or not anchor_stops_upstream and cblock is other_stop_cblock,
                                'affects_downstream_stop': anchor_stops_upstream and cblock is other_stop_cblock or not anchor_stops_upstream and cblock is anchor_stop_cblock
                            })
                fig_path = f'../output/chr19/{gene_name}.png'
                if not os.path.isfile(fig_path) and len(transcripts) > 1:
                    isoplot = IsoformPlot(transcripts)
                    isoplot.draw_all_isoforms()
                    isoplot.draw_frameshifts()
                    isoplot.savefig(fig_path)
                    plt.close(isoplot.fig)
                elif os.path.isfile(fig_path) and len(transcripts) <= 1:
                    os.remove(fig_path)
        return out

    records = []
    with mp.Pool() as p:
        t = tqdm(desc='Processing genes', total=len(gene_to_transcripts), unit='gene', file=sys.stdout)
        for result in p.imap_unordered(process_gene, gene_to_transcripts.keys()):
            records.extend(result)
            t.update()
    df = pd.DataFrame.from_records(records)
    df.to_csv(df_path, sep='\t', index=False)

df = df.fillna(value='')
display(df)

# %%
pblock_groups = df.groupby(['anchor', 'other', 'pblock'])
pblocks = pd.DataFrame(index=pd.Index(data=pblock_groups.indices.keys(), name=('anchor', 'other', 'pblock')))

pblock_attrs = pblocks.index.to_frame()['pblock'].str.extract(r'\w\((\d+):(\d+)\|(\d+):(\d+)\)').astype(int)

pblocks['category'] = pblocks.index.to_frame()['pblock'].str.get(0)
pblocks['delta_aa'] = (pblock_attrs[3] - pblock_attrs[2]) - (pblock_attrs[1] - pblock_attrs[0])

cols = ('affects_upstream_start', 'affects_downstream_start', 'affects_upstream_stop', 'affects_downstream_stop')
cblock_cat = pd.CategoricalDtype(['d', 'i', 'u', 't', 'a', 'b', '-'], ordered=True)
for col in cols:
    pblocks[col[8:]] = df['cblock'][pblock_groups[col].idxmax()].str.get(0).where(pblock_groups[col].any().array, other='-').astype(cblock_cat).array

# %%
nterm_cat = pd.CategoricalDtype(list(NTerminalChange), ordered=True)
cterm_cat = pd.CategoricalDtype(list(CTerminalChange), ordered=True)

def classify_nterm(upcat, downcat):
    if downcat in set('ut'):
        return NTerminalChange.ALTERNATIVE_ORF
    if upcat in set('di'):
        return NTerminalChange.MUTUALLY_EXCLUSIVE if downcat in set('di') else NTerminalChange.DOWNSTREAM_SHARED
    elif upcat in set('ut'):
        return NTerminalChange.UPSTREAM_SHARED if downcat in set('di') else NTerminalChange.MUTUALLY_SHARED
    else:
        return None if upcat == '-' else NTerminalChange.UNKNOWN

def classify_cterm(upcat, downcat):
    if upcat in set('di'):
        return CTerminalChange.SPLICING
    elif upcat in set('ab'):
        return CTerminalChange.FRAMESHIFT
    elif upcat in set('ut'):
        return CTerminalChange.ALTERNATIVE_ORF
    else:
        return None if downcat == '-' else CTerminalChange.UNKNOWN

pblocks['nterm'] = list(starmap(classify_nterm, zip(pblocks['upstream_start'], pblocks['downstream_start'])))
pblocks['cterm'] = list(starmap(classify_cterm, zip(pblocks['upstream_stop'], pblocks['downstream_stop'])))
pblocks['nterm'] = pblocks['nterm'].astype(nterm_cat)
pblocks['cterm'] = pblocks['cterm'].astype(cterm_cat)

# %%
pblocks['cblocks'] = pblock_groups['cblock'].apply(tuple)
pblocks['tblocks'] = pblock_groups['tblock'].unique().apply(lambda x: tuple(filter(None, x)))
pblocks['tblock_events'] = pblock_groups['events'].unique().apply(lambda x: tuple(filter(None, x)))
pblocks['events'] = pblocks['tblock_events'].apply(lambda x: set(chain.from_iterable(x)))

pblocks['complex_splicing'] = pblock_groups['complex_splicing'].agg(any)
pblocks['complex_effect'] = (pblocks['category'] == 'S') & (pblocks['cblocks'].apply(lambda cblocks: len({cblock[0] for cblock in cblocks if cblock[0] not in 'ex'})) > 1)
pblocks['frameshift'] = pblock_groups['cblock'].apply(lambda cblocks: any(cblock[0] in 'ab' for cblock in cblocks))

# %%
display(pblocks)

# %%
nterm_pblocks = pblocks[~pblocks['nterm'].isna()]
display(pd.crosstab(nterm_pblocks['upstream_start'], nterm_pblocks['downstream_start'], margins=True))
display(nterm_pblocks.value_counts('nterm'))

# %%
cterm_pblocks = pblocks[~pblocks['cterm'].isna()]
display(pd.crosstab(cterm_pblocks['upstream_stop'], cterm_pblocks['downstream_stop'], margins=True))
display(cterm_pblocks.value_counts('cterm'))

# %%
internal_pblocks = pblocks[pblocks['nterm'].isna() & pblocks['cterm'].isna()]
display(
    pd.crosstab(
        index = pblocks['complex_splicing'],
        columns = [~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()],
        margins = True
    )
)

display(
    pd.crosstab(
        index = pblocks['complex_effect'],
        columns = [~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()],
        margins = True
    )
)

display(pd.crosstab(index=pblocks['complex_splicing'], columns=pblocks['complex_effect'], margins=True))
# display(pd.crosstab(index=pblocks['complex_splicing'], columns=[pblocks['complex_effect'], ~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()], margins=True))

# %%
display(
    pd.crosstab(
        index = pblocks['nterm'],
        columns = pblocks['events'].apply(lambda x: x.intersection('Bb')).astype(bool).rename('altTSS'),
        margins = True
    )
)

# %%
