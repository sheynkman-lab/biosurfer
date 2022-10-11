# %%
import multiprocessing as mp
import sys
from itertools import chain, starmap
from operator import attrgetter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from biosurfer.core.alignments import (CodonAlignment, ProteinAlignment,
                                       SeqAlignCat, TranscriptAlignment)
from biosurfer.core.constants import APPRIS, CTerminalChange, NTerminalChange
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (ORF, GencodeTranscript, Gene,
                                                Protein,
                                                Transcript)
from biosurfer.core.splice_events import (BasicTranscriptEvent,
                                          CompoundTranscriptEvent, SpliceEvent,
                                          get_event_code)
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from matplotlib.patches import Patch, PathPatch
from more_itertools import first, one, only
from scipy.sparse import coo_matrix
from sqlalchemy import func, select
from tqdm import tqdm
import os
import pickle


def db_check(db_name):
    """Sanity check for database name.
    Args:
        db_name: User input database name

    Returns:
      Nothing
    """
    path = os.getcwd() 
    db_path = os.path.join(path, "databases")
    db_path_list = os.listdir(db_path)

    if (db_name + '.sqlite3') not in db_path_list:
        raise ValueError('\n Database ' + db_name + ' does not exist. Please create a new database or use an existing one ... \n')

   
def run_hybrid_alignment(db_name, output_path, summary):
    """ Runs hybrid alignment on input database.
    Args:
        db_name: User input database name
        output_path: Directory to write output to.

    Returns:
      Nothing
    """

    # plt.rcParams['svg.fonttype'] = 'none'
    # sns.set_style('whitegrid')

    db = Database(db_name)

    output_dir = Path(output_path).resolve()
    (output_dir/'cblock-tables').mkdir(exist_ok=True)

    # %%
    gene_to_gc_transcripts: dict[str, list[str]] = dict()

    def process_gene(gene_name: str):
        out = []
        with db.get_session() as session:
            gc_transcripts: list['GencodeTranscript'] = list(GencodeTranscript.from_names(session, gene_to_gc_transcripts[gene_name]).values())
            gc_transcripts.sort(key=attrgetter('appris'), reverse=True)
            principal = first(gc_transcripts, None)
            if not principal or not principal.protein or principal.appris is not APPRIS.PRINCIPAL or not principal.sequence:
                return out
            principal_length = principal.protein.length

            alt_transcripts: list['Transcript'] = gc_transcripts[1:]

            anchor_start_codon = principal.get_genome_coord_from_transcript_coord(principal.primary_orf.transcript_start - 1)
            anchor_stop_codon = principal.get_genome_coord_from_transcript_coord(principal.primary_orf.transcript_stop - 1)
            for alternative in alt_transcripts:
                pblocks = ()
                with ExceptionLogger(info=f'{principal}, {alternative}', callback=lambda x, y, z: alt_transcripts.remove(alternative)):
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
                        row = {
                            'anchor': principal.name,
                            'other': alternative.name,
                            'anchor_gene_id': principal.gene_id,
                            'other_gene_id': alternative.gene_id,
                            'anchor_transcript_id': principal.accession,
                            'other_transcript_id': alternative.accession,
                            'anchor_protein_id': principal.protein.accession,
                            'other_protein_id': alternative.protein.accession,
                            'pblock': str(pblock),
                            'cblock': str(cblock),
                            'cblock category': cblock.category,
                            'cblock_anchor_start': cblock.anchor_range.start,
                            'cblock_anchor_stop': cblock.anchor_range.stop,
                            'cblock_other_start': cblock.other_range.start,
                            'cblock_other_stop': cblock.other_range.stop,
                            'tblock': str(tblock),
                            'events': get_event_code(events),
                            'compound_splicing': any(set(events).intersection(compound_event.members) for compound_event in tx_aln.events if isinstance(compound_event, SpliceEvent) and len(compound_event.members) > 1),
                            'anchor_aa': principal.protein.sequence[cblock.anchor_range.start:cblock.anchor_range.stop],
                            'other_aa': alternative.protein.sequence[cblock.other_range.start:cblock.other_range.stop],
                            'affects_up_start': anchor_starts_upstream and cblock is anchor_start_cblock or not anchor_starts_upstream and cblock is other_start_cblock,
                            'affects_down_start': anchor_starts_upstream and cblock is other_start_cblock or not anchor_starts_upstream and cblock is anchor_start_cblock,
                            'affects_up_stop': anchor_stops_upstream and cblock is anchor_stop_cblock or not anchor_stops_upstream and cblock is other_stop_cblock,
                            'affects_down_stop': anchor_stops_upstream and cblock is other_stop_cblock or not anchor_stops_upstream and cblock is anchor_stop_cblock,
                            'up_start_events': '',
                            'down_start_events': '',
                            'up_stop_events': '',
                            'down_stop_events': '',
                        }
                        if cblock is anchor_start_cblock:
                            start_events = get_event_code(i.data for i in tx_aln.anchor_events.overlap(principal.primary_orf.transcript_start - 1, principal.primary_orf.transcript_start + 2) if isinstance(i.data, BasicTranscriptEvent))
                            if anchor_starts_upstream:
                                row['up_start_events'] = start_events
                            else:
                                row['down_start_events'] = start_events
                        elif cblock is other_start_cblock:
                            start_events = get_event_code(i.data for i in tx_aln.other_events.overlap(alternative.primary_orf.transcript_start - 1, alternative.primary_orf.transcript_start + 2) if isinstance(i.data, BasicTranscriptEvent))
                            if anchor_starts_upstream:
                                row['down_start_events'] = start_events
                            else:
                                row['up_start_events'] = start_events
                        if cblock is anchor_stop_cblock:
                            stop_events = get_event_code(i.data for i in tx_aln.anchor_events.overlap(principal.primary_orf.transcript_stop - 3, principal.primary_orf.transcript_stop) if isinstance(i.data, BasicTranscriptEvent))
                            if anchor_stops_upstream:
                                row['up_stop_events'] = stop_events
                            else:
                                row['down_stop_events'] = stop_events
                        elif cblock is other_stop_cblock:
                            stop_events = get_event_code(i.data for i in tx_aln.other_events.overlap(alternative.primary_orf.transcript_stop - 3, alternative.primary_orf.transcript_stop) if isinstance(i.data, BasicTranscriptEvent))
                            if anchor_stops_upstream:
                                row['down_stop_events'] = stop_events
                            else:
                                row['up_stop_events'] = stop_events
                        out.append(row)
        return out

    def process_chr(chr: str):
        """
        Read fasta file and output sequences in a list.
        Parameters
        ----------
        name : str
          Name of file containing sequences in fasta format.
        Returns
        -------
        header_lst : list
            List of headers (str)
        sequence_lst : list
            List of sequences (str)
        """
        print(f'Loading gene and transcript names for {chr}...')
        gene_to_gc_transcripts.clear()
        with db.get_session() as session:
            rows = session.execute(
                    select(Gene.name, Transcript.name, Transcript.accession, Transcript.type).
                    select_from(Protein).
                    join(Protein.orf).
                    join(ORF.transcript).
                    join(Transcript.gene).
                    where(Gene.chromosome_id == chr)
            ).all()
            for gene_name, tx_name, tx_acc, tx_type in rows:
                gc_txs = gene_to_gc_transcripts.setdefault(gene_name, [])
                if tx_type == 'gencodetranscript':
                    gc_txs.append(tx_name)

        df_path = output_dir/'cblock-tables'/f'alignment-analysis-{chr}.tsv'
        try:
            df = pd.read_csv(df_path, sep='\t')
        except:
            records = []
            # Pool()
            t = tqdm(desc='Processing genes', total=len(gene_to_gc_transcripts), unit='gene', file=sys.stdout)
            for result in map(process_gene, gene_to_gc_transcripts.keys()):
                records.extend(result)
                t.update()
            df = pd.DataFrame.from_records(records)
            df.to_csv(df_path, sep='\t', index=False)
        return df

    chrs = [f'chr{i}' for i in list(range(1, 23)) + ['X']]
    df = pd.concat(
        (process_chr(chr) for chr in chrs),
        keys = chrs,
        names = ['chr', 'row']
    ).fillna(value='').reset_index().drop(columns='row')

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

    # %%
    pblock_groups = df.groupby(['chr', 'anchor', 'other', 'pblock'])
    pblock_stat = df.groupby(['pblock'])['pblock'].count()

    pblocks = pd.DataFrame(index=pd.Index(data=pblock_groups.indices.keys(), name=('chr', 'anchor', 'other', 'pblock')))

    pblock_attrs = pblocks.index.to_frame()['pblock'].str.extract(r'\w\((\d+):(\d+)\|(\d+):(\d+)\)').astype(int)

    pblock_cat = pd.CategoricalDtype(['I', 'D', 'S'], ordered=True)

    pblocks['category'] = pblocks.index.to_frame()['pblock'].str.get(0).astype(pblock_cat)

    pblocks['length_change'] = (pblock_attrs[3] - pblock_attrs[2]) - (pblock_attrs[1] - pblock_attrs[0])
    pblocks['anchor_length'] = [protein_lengths[anchor] for anchor in pblocks.reset_index()['anchor']]
    pblocks['anchor_relative_length_change'] = pblocks['length_change'] / pblocks['anchor_length']

    # %%
    cblock_cat = pd.CategoricalDtype(['d', 'i', 'u', 't', 'a', 'b', '-'], ordered=True)
    for col in ('affects_up_start', 'affects_down_start', 'affects_up_stop', 'affects_down_stop'):
        indices = pblock_groups[col].idxmax()
        pblocks[col[8:] + '_cblock'] = df['cblock'][indices].str.get(0).where(pblock_groups[col].any().array, other='-').astype(cblock_cat).array
        pblocks[col[8:] + '_cblock_events'] = df['events'][indices].where(pblock_groups[col].any().array, other='').array

    for col in ('up_start_events', 'down_start_events', 'up_stop_events', 'down_stop_events'):
        pblocks[col] = pblock_groups[col].max()

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

    pblocks['nterm'] = list(starmap(classify_nterm, zip(pblocks['up_start_cblock'], pblocks['down_start_cblock'])))
    pblocks['cterm'] = list(starmap(classify_cterm, zip(pblocks['up_stop_cblock'], pblocks['down_stop_cblock'])))
    pblocks['nterm'] = pblocks['nterm'].astype(nterm_cat)
    pblocks['cterm'] = pblocks['cterm'].astype(cterm_cat)
    pblocks['internal'] = pblocks['nterm'].isna() & pblocks['cterm'].isna()

    # %%
    pblocks['cblocks'] = pblock_groups['cblock'].apply(tuple)
    pblocks['tblocks'] = pblock_groups['tblock'].unique().apply(lambda x: tuple(filter(None, x)))
    pblocks['tblock_events'] = pblock_groups['events'].unique().apply(lambda x: tuple(filter(None, x)))
    pblocks['events'] = pblocks['tblock_events'].apply(lambda x: frozenset(chain.from_iterable(x)))

    pblocks['compound_splicing'] = pblock_groups['compound_splicing'].agg(any)
    pblocks['frameshift'] = pblock_groups['cblock'].apply(lambda cblocks: any(cblock[0] in 'ab' for cblock in cblocks))
    pblocks['split_ends'] = pblock_groups['cblock'].apply(lambda cblocks: any(cblock[0] in 'ex' for cblock in cblocks))


    print("\n ----- Hybrid alignment done .... \n")
    print("------------------------------------ \n")
    # %%
    if summary:

        print("\n ----- Printing hybrid alignment statistic .... \n")
        print("------------------------------------ \n")
        #%%
        newpath = str(output_dir) + '/plots/'
        if not os.path.exists(newpath):
            os.makedirs(newpath)

        display(df)
        display(pblocks)

        #%% 
        # Plotting number of alt isoforms for N-termini
        nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != NTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.ALTERNATIVE_ORF)].copy()
        nterm_pblocks['nterm'] = nterm_pblocks['nterm'].cat.remove_unused_categories()
        nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: x.intersection('BbPp')).astype(bool)

        print("\n ------------------------------------ \n")
        display(pd.crosstab(nterm_pblocks['up_start_cblock'], nterm_pblocks['down_start_cblock'], margins=True))
        print("\n ------------------------------------ \n")

        nterm_palette = dict(zip(NTerminalChange, sns.color_palette('viridis_r', n_colors=5)[:-1]))

        nterm_fig = plt.figure(figsize=(3, 5))
        ax = sns.countplot(
            data = nterm_pblocks,
            y = 'nterm',
            order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED, NTerminalChange.UPSTREAM_SHARED, NTerminalChange.MUTUALLY_SHARED),
            palette = nterm_palette,
            linewidth = 0,
            saturation = 1,
        )
        ax.set(xlabel='Number of alternative isoforms', ylabel='Alternate N-termini classifications', yticklabels=['MXS', 'SDS', 'SUS', 'MSS']) 
        plt.savefig(str(output_dir) + '/plots/nterm-class-counts.svg', dpi=200, facecolor=None)

        #%%
        # Plotting chnage in N-terminal length
        nterm_length_order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED, NTerminalChange.MUTUALLY_SHARED)
        nterm_length_fig = plt.figure(figsize=(6, 4))
        ax = sns.violinplot(
            data = nterm_pblocks,
            x = 'anchor_relative_length_change',
            y = 'nterm',
            order = nterm_length_order,
            gridsize = 200,
            palette = nterm_palette,
            saturation = 1,
            scale = 'area',
        )

        xmax = max(ax.get_xlim())
        ymin, ymax = ax.get_ylim()
        ax.vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linewidth=1, linestyle=':')
        ax.set(xlim=(-1, 1), ylim=(ymin, ymax), xlabel='Change in N-terminal length (fraction of anchor isoform length)', ylabel=None, yticklabels=['MXS', 'SDS', 'MSS'])
        plt.savefig(str(output_dir) + '/plots/nterm-length-change-dist.svg', dpi=200, facecolor=None)
        #%%
        # Number of N-termini alternative isoforms
        tss_fig = plt.figure(figsize=(3, 10)) # Added height
        ax = sns.countplot(
            data = nterm_pblocks,
            x = 'nterm',
            palette = nterm_palette,
            saturation = 1,
            order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED),
        )
        sns.countplot(
            ax = ax,
            data = nterm_pblocks[nterm_pblocks['altTSS']],
            x = 'nterm',
            order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED),
            fill = False,
            edgecolor = 'w',
            hatch = '//',
        )
        ax.legend(loc=(0, 1), frameon=False, handles=[Patch(facecolor='k', edgecolor='w', hatch='///'), Patch(facecolor='k')], labels=['driven by alternate TSS', 'driven by 5\' UTR splicing'])
        ax.set(ylabel='Number of alternative isoforms', xlabel=None, xticklabels=['MXS', 'SDS'])
        plt.savefig(str(output_dir) + '/plots/nterm-altTSS-counts.svg', dpi=200, facecolor=None)

        #%%
        # Number of C-termini alternative isoforms
        cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'] != NTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.UNKNOWN)].copy()
        cterm_pblocks['cterm'] = cterm_pblocks['cterm'].cat.remove_unused_categories()
        cterm_pblocks['APA'] = cterm_pblocks['events'].apply(lambda x: x.intersection('BbPp')).astype(bool)

        print("\n ------------------------------------ \n")
        display(pd.crosstab(cterm_pblocks['up_stop_cblock'], cterm_pblocks['down_stop_cblock'], margins=True))
        print("\n ------------------------------------ \n")

        cterm_splice_palette = sns.color_palette('RdPu_r', n_colors=3)
        cterm_frameshift_palette = sns.color_palette('YlOrRd_r', n_colors=4)
        cterm_palette = [cterm_splice_palette[0], cterm_frameshift_palette[0]]

        cterm_fig = plt.figure(figsize=(3, 2))
        ax = sns.countplot(
            data = cterm_pblocks,
            y = 'cterm',
            order = (CTerminalChange.SPLICING, CTerminalChange.FRAMESHIFT),
            palette = cterm_palette,
            saturation = 1,
            linewidth = 0,
        )
        ax.set(xlabel='Number of alternative isoforms', ylabel=None, yticklabels=[])
        plt.savefig(str(output_dir) + '/plots/cterm-class-counts.svg', dpi=200, facecolor=None)

        #%% 
        # Plotting all the C-termini sub categories
        cterm_pblock_events = cterm_pblocks['up_stop_events'].combine(cterm_pblocks['down_stop_events'], lambda x, y: (x, y))
        single_ATE = (cterm_pblocks['cterm'] == CTerminalChange.SPLICING) & cterm_pblocks['tblock_events'].isin({('B', 'b'), ('b', 'B')})
        cterm_splice_subcats = pd.DataFrame(
            {
                'exon extension introduces termination': cterm_pblocks['up_stop_events'].isin({'P', 'I', 'D'}),
                'alternative terminal exon(s)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}),
                'poison exon': cterm_pblocks['up_stop_events'] == 'E',
                'other': [True for _ in cterm_pblocks.index]
            }
        )
        cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

        cterm_splice_palette_dict = dict(zip(
            cterm_splice_subcats.columns,
            cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + ['#bbbbbb']
        ))

        splice_subcat_order = tuple(cterm_splice_subcats.keys())

        cterm_splice_fig, axs = plt.subplots(1, 2, figsize=(10, 6))
        sns.countplot(
            ax = axs[0],
            data = cterm_pblocks[cterm_pblocks['cterm'] == CTerminalChange.SPLICING],
            y = 'splice_subcat',
            order = splice_subcat_order,
            palette = cterm_splice_palette_dict,
            saturation = 1,
            linewidth = 0,
        )
        axs[0].set(xlabel='Number of alternative isoforms', ylabel=None)

        sns.violinplot(
            ax = axs[1],
            data = cterm_pblocks[cterm_pblocks['cterm'] == CTerminalChange.SPLICING],
            x = 'anchor_relative_length_change',
            y = 'splice_subcat',
            order = splice_subcat_order,
            palette = cterm_splice_palette_dict,
            saturation = 1,
            gridsize = 200,
            scale = 'area',
        )
        xmax = max(axs[1].get_xlim())
        ymin, ymax = axs[1].get_ylim()
        axs[1].vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linewidth=1, linestyle=':')
        axs[1].set(xlim=(-1, 1), ylim=(ymin, ymax), xlabel='change in C-terminal length (fraction of anchor isoform length)', ylabel=None, yticklabels=[])

        plt.savefig(str(output_dir) + '/plots/cterm-splicing-subcats.svg', dpi=200, facecolor=None, bbox_inches='tight')

        #%%
        # Plotting all the C-termini frameshift subcategories
        cterm_frame_subcats = pd.DataFrame(
        {
            'exon': cterm_pblocks['up_stop_cblock_events'].isin({'E', 'e'}),
            'acceptor': cterm_pblocks['up_stop_cblock_events'].isin({'A', 'a'}),
            'donor': cterm_pblocks['up_stop_cblock_events'].isin({'D', 'd'}),
            'intron': cterm_pblocks['up_stop_cblock_events'].isin({'I', 'i'}),
            'other': [True for _ in cterm_pblocks.index]
        }
        )
        cterm_pblocks['frame_subcat'] = cterm_frame_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_frame_subcats.columns, ordered=True))

        cterm_frameshift_palette_dict = dict(zip(
            cterm_frame_subcats.columns,
            cterm_frameshift_palette + ['#bbbbbb']
        ))

        frame_subcat_order = cterm_pblocks[cterm_pblocks['cterm'] == CTerminalChange.FRAMESHIFT]['frame_subcat'].value_counts().index

        cterm_frameshift_fig, axs = plt.subplots(1, 2, figsize=(10, 6))
        sns.countplot(
            ax = axs[0],
            data = cterm_pblocks[cterm_pblocks['cterm'] == CTerminalChange.FRAMESHIFT],
            y = 'frame_subcat',
            order = frame_subcat_order,
            palette = cterm_frameshift_palette_dict,
            saturation = 1,
            linewidth = 0,
        )
        axs[0].set(xlabel='Number of alternative isoforms', ylabel=None)

        sns.violinplot(
            ax = axs[1],
            data = cterm_pblocks[cterm_pblocks['cterm'] == CTerminalChange.FRAMESHIFT],
            x = 'anchor_relative_length_change',
            y = 'frame_subcat',
            order = frame_subcat_order,
            palette = cterm_frameshift_palette_dict,
            saturation = 1,
            scale = 'area'
        )
        xmax = max(axs[1].get_xlim())
        ymin, ymax = axs[1].get_ylim()
        axs[1].vlines(x=0, ymin=ymin, ymax=ymax, color='#444444', linewidth=1, linestyle=':')
        axs[1].set(xlim=(-1, 1), ylim=(ymin, ymax), xlabel='change in C-terminal length (fraction of anchor isoform length)', ylabel=None, yticklabels=[])

        plt.savefig(str(output_dir) + '/plots/cterm-frameshift-subcats.svg', dpi=200, facecolor=None, bbox_inches='tight')

        #%%
        # Plotting internal regions pblock counts for alt isoforms
        cterm_event_counts = cterm_pblocks.groupby('cterm').events.value_counts().rename('count')
        cterm_examples = cterm_pblocks.groupby(['cterm', 'events']).sample(1, random_state=329).join(cterm_event_counts, on=['cterm', 'events']).sort_values(['cterm', 'count'], ascending=False).set_index(['cterm', 'events'])

        internal_pblocks = (
            pblocks[pblocks['nterm'].isna() & pblocks['cterm'].isna()].
            drop(columns=[col for col in pblocks.columns if 'start' in col or 'stop' in col]).
            copy()
        )
        internal_pblocks['category'] = (
            internal_pblocks['cblocks'].
            apply(lambda cblocks: ''.join(cblock[0] for cblock in cblocks)).
            str.replace(r'[ex]', '', regex=True).
            map({'d': 'D', 'i': 'I'}).
            fillna('S')
        )
        internal_pblock_counts = internal_pblocks.reset_index(level=3).groupby(['anchor', 'other']).agg(pblocks=('pblock', 'count'))

        internal_pblock_counts_fig = plt.figure(figsize=(6, 4))
        ax = sns.countplot(data=internal_pblock_counts, x='pblocks', palette='Blues_r')
        ax.set(xlabel='Number of alternative internal regions', ylabel='Number of alternative isoforms')
        plt.savefig(str(output_dir) + '/plots/internal-pblock-counts.svg', dpi=200, facecolor=None)

        # %%
        # Plotting internal pblock AS events for interal regions
        internal_cat_palette = {'D': '#f800c0', 'I': '#00c0f8', 'S': '#f8c000'}
        internal_event_palette = {
            'Intron': '#e69138',
            'Alt. donor': '#6aa84f',
            'Alt. acceptor': '#8a4ea7',
            'Single exon': '#3d85c6',
            'Mutually exclusive exons': '#255179',
            'Compound': '#888888'
        }

        internal_subcats = pd.DataFrame(
            {
                'Intron': internal_pblocks['tblock_events'].isin({('I',), ('i',)}),
                'Alt. donor': internal_pblocks['tblock_events'].isin({('D',), ('d',)}),
                'Alt. acceptor': internal_pblocks['tblock_events'].isin({('A',), ('a',)}),
                'Single exon': internal_pblocks['tblock_events'].isin({('E',), ('e',)}),
                'Mutually exclusive exons': internal_pblocks['tblock_events'].isin({('E', 'e'), ('e', 'E')}),
                'Compound': [True for _ in internal_pblocks.index]
            }
        )
        internal_pblocks['splice event'] = internal_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(internal_subcats.columns, ordered=True))

        internal_pblocks_fig = plt.figure(figsize=(4, 4))
        ax = sns.countplot(
            data = internal_pblocks.sort_values('category', ascending=True),
            y = 'splice event',
            hue = 'category',
            palette = internal_cat_palette,
            saturation = 1,
            dodge = True,
        )
        plt.legend(loc='lower right', labels=['Deletion', 'Insertion', 'Substitution'])
        ax.set(xlabel='Number of alternative internal regions', ylabel=None)
        plt.savefig(str(output_dir) + '/plots/internal-pblock-events.svg', dpi=200, facecolor=None, bbox_inches='tight')

        #%%
        # Plotting splice events for internal pblocks
        internal_pblocks_split_fig = plt.figure(figsize=(4, 4))
        ax = sns.countplot(
            data = internal_pblocks.sort_values('category', ascending=True),
            y = 'splice event',
            palette = internal_event_palette,
            saturation = 1,
        )
        sns.countplot(
            ax = ax,
            data = internal_pblocks[internal_pblocks.split_ends].sort_values('category', ascending=True),
            y = 'splice event',
            fill = False,
            edgecolor = 'w',
            hatch = '///',
        )
        plt.legend(loc='lower right', handles=[Patch(facecolor='k', edgecolor='w', hatch='///'), Patch(facecolor='k', edgecolor='w')], labels=['split codons', 'no split codons'])
        ax.set(xlabel='Number of alternative internal regions', ylabel=None)
        plt.savefig(str(output_dir) + '/plots/internal-pblock-events-split.svg', dpi=200, facecolor=None, bbox_inches='tight')

        #%%
        # Plotting compound splicing events for pblocks in internal region
        internal_compound_pblocks = internal_pblocks[internal_pblocks['splice event'] == 'Compound'].copy()
        internal_compound_subcats = pd.DataFrame(
            {
                'Multi-exon skip': internal_compound_pblocks['events'] == {'e'},
                'Exon skipping + alt donor/acceptor': internal_compound_pblocks['events'].isin({
                    frozenset('de'),
                    frozenset('De'),
                    frozenset('ea'),
                    frozenset('eA'),
                    frozenset('dea'),
                    frozenset('Dea'),
                    frozenset('deA'),
                    frozenset('DeA'),
                }),
                'Other': [True for _ in internal_compound_pblocks.index]
            }
        )
        internal_compound_pblocks['compound_subcat'] = internal_compound_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(internal_compound_subcats.columns, ordered=True))

        internal_pblocks_compound_fig = plt.figure(figsize=(6, 1.5))
        ax = sns.countplot(
                data = internal_compound_pblocks,
                y = 'compound_subcat',
                palette = 'Blues_r',
                saturation = 1,
                linewidth = 0,
        )
        ax.set(xlabel='Number of alternative internal regions', ylabel=None)
        plt.savefig(str(output_dir) + '/plots/internal-pblock-compound-events.svg', dpi=200, facecolor=None, bbox_inches='tight')


    # %%

    # Pickle the pblocks datafrom for vis.
    # pickle_out = open("biosurfer.pickle","wb")
    # pickle.dump(pblocks, pickle_out)
    # pickle.dump(output_path, pickle_out)
    # pickle_out.close()
