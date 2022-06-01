# %%
from functools import reduce
from math import ceil
import multiprocessing as mp
import os
from pathlib import Path
import sys
from itertools import chain, groupby, starmap
from operator import attrgetter, or_

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
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
from matplotlib.patches import Patch
from more_itertools import first, one
from sqlalchemy import select, func
from tqdm import tqdm

db = Database('gencode')
output_dir = Path('../output/alignment-analysis')

# %%
gene_to_transcripts: dict[str, list[str]] = dict()

def process_gene(gene_name: str):
    out = []
    with db.get_session() as session:
        transcripts: list['GencodeTranscript'] = list(GencodeTranscript.from_names(session, gene_to_transcripts[gene_name]).values())
        transcripts.sort(key=attrgetter('appris'), reverse=True)
        principal = transcripts[0]
        principal_length = principal.protein.length
        if principal.appris is APPRIS.PRINCIPAL and principal.sequence:
            anchor_start_codon = principal.get_genome_coord_from_transcript_coord(principal.primary_orf.transcript_start - 1)
            anchor_stop_codon = principal.get_genome_coord_from_transcript_coord(principal.primary_orf.transcript_stop - 1)
            for alternative in transcripts[1:]:
                pblocks = ()
                with ExceptionLogger(info=f'{principal}, {alternative}', callback=lambda x, y, z: transcripts.remove(alternative)):
                    other_start_codon = alternative.get_genome_coord_from_transcript_coord(alternative.primary_orf.transcript_start - 1)
                    other_stop_codon = alternative.get_genome_coord_from_transcript_coord(alternative.primary_orf.transcript_stop - 1)
                    anchor_starts_upstream = anchor_start_codon <= other_start_codon
                    anchor_stops_upstream = anchor_stop_codon <= other_stop_codon

                    alternative_length = alternative.protein.length
                    tx_aln = TranscriptAlignment.from_transcripts(principal, alternative)
                    cd_aln = CodonAlignment.from_proteins(principal.protein, alternative.protein)
                    pr_aln = ProteinAlignment.from_proteins(principal.protein, alternative.protein)
                    pblocks = pr_aln.blocks
                    anchor_start_cblock = one(cd_aln.anchor_blocks.at(0)).data
                    other_start_cblock = one(cd_aln.other_blocks.at(0)).data
                    anchor_stop_cblock = one(cd_aln.anchor_blocks.at(principal_length - 1)).data
                    other_stop_cblock = one(cd_aln.other_blocks.at(alternative_length - 1)).data
                for pblock in pblocks:
                    if pblock.category is SeqAlignCat.MATCH:
                        continue
                    for cblock in pr_aln.pblock_to_cblocks[pblock]:
                        tblock = cd_aln.cblock_to_tblock.get(cblock)
                        events = tx_aln.block_to_events.get(tblock, ())
                        out.append({
                            'anchor': principal.name,
                            'other': alternative.name,
                            'pblock': str(pblock),
                            'cblock': str(cblock),
                            'tblock': str(tblock),
                            'events': get_event_code(events),
                            'compound_splicing': any(set(events).intersection(compound_event.members) for compound_event in tx_aln.events if isinstance(compound_event, SpliceEvent) and len(compound_event.members) > 1),
                            'affects_upstream_start': anchor_starts_upstream and cblock is anchor_start_cblock or not anchor_starts_upstream and cblock is other_start_cblock,
                            'affects_downstream_start': anchor_starts_upstream and cblock is other_start_cblock or not anchor_starts_upstream and cblock is anchor_start_cblock,
                            'affects_upstream_stop': anchor_stops_upstream and cblock is anchor_stop_cblock or not anchor_stops_upstream and cblock is other_stop_cblock,
                            'affects_downstream_stop': anchor_stops_upstream and cblock is other_stop_cblock or not anchor_stops_upstream and cblock is anchor_stop_cblock
                        })
            # fig_path = output_dir/f'chr19/{gene_name}.png'
            # if not fig_path.isfile() and len(transcripts) > 1:
            #     isoplot = IsoformPlot(transcripts)
            #     isoplot.draw_all_isoforms()
            #     isoplot.draw_frameshifts()
            #     isoplot.savefig(fig_path)
            #     plt.close(isoplot.fig)
            # elif fig_path.isfile() and len(transcripts) <= 1:
            #     os.remove(fig_path)
    return out

def process_chr(chr: str):
    print(f'Loading gene and transcript names for {chr}...')
    gene_to_transcripts.clear()
    with db.get_session() as session:
        rows = session.execute(
                select(Gene.name, Transcript.name).
                select_from(Protein).
                join(Protein.orf).
                join(ORF.transcript).
                join(Transcript.gene).
                where(Gene.chromosome_id == chr)
        ).all()
        for gene_name, tx_name in rows:
            gene_to_transcripts.setdefault(gene_name, []).append(tx_name)

    df_path = output_dir/f'alignment-analysis-{chr}.tsv'
    try:
        df = pd.read_csv(df_path, sep='\t')
    except:
        records = []
        with mp.Pool() as p:
            t = tqdm(desc='Processing genes', total=len(gene_to_transcripts), unit='gene', file=sys.stdout)
            for result in p.imap_unordered(process_gene, gene_to_transcripts.keys()):
                records.extend(result)
                t.update()
        df = pd.DataFrame.from_records(records)
        df.to_csv(df_path, sep='\t', index=False)
    return df

df = pd.concat((process_chr(f'chr{i}') for i in list(range(19, 20))), ignore_index=True).fillna(value='')
display(df)

with db.get_session() as session:
    protein_lengths = {
        tx_name: protein_length
        for tx_name, protein_length in session.execute(
            select(Transcript.name, func.length(Protein.sequence)).
            join(Protein.orf).
            join(ORF.transcript).
            where(Transcript.name.in_(list(chain(df['anchor'].unique(), df['other'].unique()))))
        ).all()
    }

# sys.exit()

# %%
pblock_groups = df.groupby(['anchor', 'other', 'pblock'])
pblocks = pd.DataFrame(index=pd.Index(data=pblock_groups.indices.keys(), name=('anchor', 'other', 'pblock')))

pblock_attrs = pblocks.index.to_frame()['pblock'].str.extract(r'\w\((\d+):(\d+)\|(\d+):(\d+)\)').astype(int)

pblocks['category'] = pblocks.index.to_frame()['pblock'].str.get(0)
pblocks['length_change'] = (pblock_attrs[3] - pblock_attrs[2]) - (pblock_attrs[1] - pblock_attrs[0])
pblocks['anchor_length'] = [protein_lengths[anchor] for anchor in pblocks.reset_index()['anchor']]
pblocks['anchor_relative_length_change'] = pblocks['length_change'] / pblocks['anchor_length']

# %%
cols = ('affects_upstream_start', 'affects_downstream_start', 'affects_upstream_stop', 'affects_downstream_stop')
cblock_cat = pd.CategoricalDtype(['d', 'i', 'u', 't', 'a', 'b', '-'], ordered=True)
for col in cols:
    pblocks[col[8:]] = df['cblock'][pblock_groups[col].idxmax()].str.get(0).where(pblock_groups[col].any().array, other='-').astype(cblock_cat).array

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
pblocks['internal'] = pblocks['nterm'].isna() & pblocks['cterm'].isna()

# %%
pblocks['cblocks'] = pblock_groups['cblock'].apply(tuple)
pblocks['tblocks'] = pblock_groups['tblock'].unique().apply(lambda x: tuple(filter(None, x)))
pblocks['tblock_events'] = pblock_groups['events'].unique().apply(lambda x: tuple(filter(None, x)))
pblocks['events'] = pblocks['tblock_events'].apply(lambda x: frozenset(chain.from_iterable(x)))

pblocks['compound_splicing'] = pblock_groups['compound_splicing'].agg(any)
# pblocks['complex_effect'] = (pblocks['category'] == 'S') & (pblocks['cblocks'].apply(lambda cblocks: len({cblock[0] for cblock in cblocks if cblock[0] not in 'ex'})) > 1)
pblocks['frameshift'] = pblock_groups['cblock'].apply(lambda cblocks: any(cblock[0] in 'ab' for cblock in cblocks))

# %%
def get_genes(pblocks_slice: 'pd.DataFrame'):
    return pblocks_slice.reset_index()['anchor'].str.extract(r'(\w+)-\d+', expand=False).unique()

with open('../output/genes-mxs.txt', 'w') as f:
    f.writelines(s+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.MUTUALLY_EXCLUSIVE]))

with open('../output/genes-sds.txt', 'w') as f:
    f.writelines(s+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.DOWNSTREAM_SHARED]))

with open('../output/genes-sus.txt', 'w') as f:
    f.writelines(s+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.UPSTREAM_SHARED]))

with open('../output/genes-mss.txt', 'w') as f:
    f.writelines(s+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.MUTUALLY_SHARED]))

# %%
display(pblocks)

# %%
nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != NTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.ALTERNATIVE_ORF)].copy()
nterm_pblocks['nterm'] = nterm_pblocks['nterm'].cat.remove_unused_categories()
nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: x.intersection('BbPp')).astype(bool)
display(pd.crosstab(nterm_pblocks['upstream_start'], nterm_pblocks['downstream_start'], margins=True))

nterm_palette = sns.color_palette('viridis_r', n_colors=4)

nterm_fig = plt.figure(figsize=(3, 4))
ax = sns.countplot(
    data = nterm_pblocks,
    y = 'nterm',
    order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED, NTerminalChange.UPSTREAM_SHARED, NTerminalChange.MUTUALLY_SHARED),
    palette = nterm_palette
)
ax.set(xlabel='# of alternative isoforms', ylabel=None, yticklabels=[])
plt.savefig(output_dir/'nterm-class-counts.png', dpi=200, facecolor=None)

# %%
nterm_length_fig = plt.figure(figsize=(6, 4))
ax = sns.stripplot(
    data = nterm_pblocks,
    x = 'anchor_relative_length_change',
    y = 'nterm',
    order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED),
    jitter = 0.25,
    size = 2,
    palette = nterm_palette[:2]
)
xmax = max(ax.get_xlim())
ymin, ymax = sorted(ax.get_ylim())
ax.vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linestyle=':')
ax.set(xlim=(-1, ceil(xmax)), xlabel='change in N-terminal length (fraction of anchor isoform length)', ylabel=None, yticklabels=['MXS', 'SDS'])
plt.savefig(output_dir/'nterm-length-change-dist.png', dpi=200, facecolor=None)

# %%
tss_fig = plt.figure(figsize=(6, 4))
ax = sns.countplot(
    data = nterm_pblocks,
    y = 'nterm',
    palette = nterm_palette,
    order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED)
)
sns.countplot(
    ax = ax,
    data = nterm_pblocks[nterm_pblocks['altTSS']],
    y = 'nterm',
    order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED),
    fill = False,
    edgecolor = 'w',
    hatch = '//',
)
ax.legend(handles=[Patch(facecolor='gray', edgecolor='w', hatch='///'), Patch(facecolor='gray')], labels=['driven by alternate TSS', 'driven by 5\' UTR splicing'])
ax.set(xlabel='# of alternative isoforms', ylabel=None, yticklabels=['MXS', 'SDS'])
plt.savefig(output_dir/'nterm-altTSS-counts.png', dpi=200, facecolor=None)

# %%
cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'] != NTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.UNKNOWN)].copy()
cterm_pblocks['cterm'] = cterm_pblocks['cterm'].cat.remove_unused_categories()
cterm_pblocks['APA'] = cterm_pblocks['altTSS'] = cterm_pblocks['events'].apply(lambda x: x.intersection('BbPp')).astype(bool)
# cterm_pblocks['ATE'] = cterm_pblocks['events'] == {'b', 'B'}
cterm_pblocks['ATE'] = cterm_pblocks.tblock_events.isin([('B', 'b'), ('b', 'B')])
cterm_pblocks['EXIT'] = cterm_pblocks['events'] == {'b', 'P'}  # TODO: this doesn't cover all kinds of EXIT, but it seems to be the most common
cterm_pblocks['exon'] = (cterm_pblocks['events'] == {'E'}) | (cterm_pblocks['events'] == {'e'})
cterm_pblocks['intron'] = (cterm_pblocks['events'] == {'I'}) | (cterm_pblocks['events'] == {'i'})
display(pd.crosstab(cterm_pblocks['upstream_stop'], cterm_pblocks['downstream_stop'], margins=True))

cterm_palette = sns.color_palette('rocket', n_colors=2)

cterm_fig = plt.figure(figsize=(3, 4))
ax = sns.countplot(
    data = cterm_pblocks,
    y = 'cterm',
    order = (CTerminalChange.SPLICING, CTerminalChange.FRAMESHIFT),
    palette = cterm_palette
)
ax.set(xlabel='# of alternative isoforms', ylabel=None, yticklabels=[])
plt.savefig(output_dir/'cterm-class-counts.png', dpi=200, facecolor=None)

# %%
cterm_length_fig = plt.figure(figsize=(6, 4))
ax = sns.stripplot(
    data = cterm_pblocks,
    x = 'anchor_relative_length_change',
    y = 'cterm',
    order = (CTerminalChange.SPLICING, CTerminalChange.FRAMESHIFT),
    jitter = 0.25,
    size = 2,
    palette = cterm_palette
)
xmax = max(ax.get_xlim())
ymin, ymax = sorted(ax.get_ylim())
ax.vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linestyle=':')
ax.set(xlim=(-1, ceil(xmax)), xlabel='change in C-terminal length (fraction of anchor isoform length)', ylabel=None, yticklabels=['splice', 'frame'])
plt.savefig(output_dir/'cterm-length-change-dist.png', dpi=200, facecolor=None)

# %%
cterm_event_fig = plt.figure(figsize=(6, 4))
ax = sns.countplot(
    data = cterm_pblocks,
    y = 'cterm',
    palette = cterm_palette,
    order = (CTerminalChange.SPLICING, CTerminalChange.FRAMESHIFT)
)
subcats = ['ATE', 'EXIT', 'exon', 'intron']
hatches = ['o.', '..', 'xx', '//']
for i in reversed(range(len(subcats))):
    sns.countplot(
        ax = ax,
        data = cterm_pblocks[reduce(or_, (cterm_pblocks[subcat] for subcat in subcats[:i+1]))],
        y = 'cterm',
        order = (CTerminalChange.SPLICING, CTerminalChange.FRAMESHIFT),
        palette = cterm_palette,
        edgecolor = 'w',
        hatch = hatches[i],
        label = subcats[i]
    )
ax.legend(ncol=2)
# ax.legend(handles=[Patch(facecolor='gray', edgecolor='w', hatch='///'), Patch(facecolor='gray')], labels=['driven by alternate TSS', 'driven by 5\' UTR splicing'])
ax.set(xlabel='# of alternative isoforms', ylabel=None, yticklabels=['splice', 'frame'])
plt.savefig(output_dir/'cterm-event-counts.png', dpi=200, facecolor=None)

# %%
cterm_event_counts = cterm_pblocks.groupby('cterm').events.value_counts().rename('count')
cterm_examples = cterm_pblocks.groupby(['cterm', 'events']).sample(1, random_state=329).join(cterm_event_counts, on=['cterm', 'events']).sort_values(['cterm', 'count'], ascending=False).set_index(['cterm', 'events'])

# %%
# internal_pblocks = pblocks[pblocks['nterm'].isna() & pblocks['cterm'].isna()]
# display(
#     pd.crosstab(
#         index = pblocks['complex_splicing'],
#         columns = [~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()],
#         margins = True
#     )
# )

# display(
#     pd.crosstab(
#         index = pblocks['complex_effect'],
#         columns = [~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()],
#         margins = True
#     )
# )

# display(pd.crosstab(index=pblocks['complex_splicing'], columns=pblocks['complex_effect'], margins=True))
# display(pd.crosstab(index=pblocks['complex_splicing'], columns=[pblocks['complex_effect'], ~pblocks['nterm'].isna(), ~pblocks['cterm'].isna()], margins=True))

# %%
