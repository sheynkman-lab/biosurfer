# %%
import multiprocessing as mp
import sys
from itertools import chain, starmap
from operator import attrgetter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from biosurfer.core.alignments import (CodonAlignment, ProteinAlignment,
                                       SeqAlignCat, TranscriptAlignment)
from biosurfer.core.constants import APPRIS, CTerminalChange, NTerminalChange
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import (ORF, GencodeTranscript, Gene, PacBioTranscript,
                                                Protein,
                                                Transcript)
from biosurfer.core.splice_events import (BasicTranscriptEvent,
                                          CompoundTranscriptEvent, SpliceEvent,
                                          get_event_code)
from IPython.display import display
from more_itertools import first, one, only
from sqlalchemy import func, select, and_
from tqdm import tqdm
import os
import pickle



def run_hybrid_alignment_for_all_genes(db_name, output_dir: 'Path', gencode: bool = False, gene_to_anchor_tx: dict[str, str] = None):
    """ Runs hybrid alignment on input database and saves tables of p-blocks and c-blocks.
    Args:
        db_name: User input database name
        output_path: Directory to write output to.

    Returns:
      Nothing
    """

    # plt.rcParams['svg.fonttype'] = 'none'
    # sns.set_style('whitegrid')

    cblock_dir = output_dir/'cblock-tables'
    log_dir = output_dir/'alignment-errors'
    for dir in (cblock_dir, log_dir):
        dir.mkdir(exist_ok=True)

    # db = Database(db_name)
    cblock_df = get_cblocks(db_name, cblock_dir, log_dir, gencode, gene_to_anchor_tx)
    display(cblock_df)
    # pblock_df = generate_pblock_table(cblock_df)


def process_chr(chr: str, db_name: str, log_file: 'Path', gencode: bool, gene_to_anchor_tx: dict[str, str]):
    db = Database(db_name)
    gene_to_gc_transcripts: dict[str, list[str]] = dict()
    gene_to_pb_transcripts: dict[str, list[str]] = dict()
    open(log_file, 'w').close()  # clear log_file contents

    def process_gene(gene_name: str):
        out = []
        with db.get_session() as session:
            gc_transcripts: dict[str, 'GencodeTranscript'] = GencodeTranscript.from_names(session, gene_to_gc_transcripts[gene_name])
            pb_transcripts: dict[str, 'PacBioTranscript'] = {tx.accession: tx for tx in PacBioTranscript.from_accessions(session, gene_to_pb_transcripts[gene_name]).values()}
            
            # choose anchor isoform for gene
            if gene_to_anchor_tx:
                anchor_id = gene_to_anchor_tx.get(gene_name, None)
                if anchor_id in gc_transcripts:
                    anchor = gc_transcripts[anchor_id]
                elif anchor_id in pb_transcripts:
                    anchor = pb_transcripts[anchor_id]
                else:
                    anchor = None
            elif gc_transcripts:  # by default, use APPRIS principal as anchor
                anchor = max(gc_transcripts.values(), key=attrgetter('appris'))
                if anchor.appris is not APPRIS.PRINCIPAL:
                    anchor = None
            else:
                anchor = None
            if not anchor or not anchor.protein or not anchor.sequence:
                return out
            
            principal_length = anchor.protein.length
            anchor_start_codon = anchor.get_genome_coord_from_transcript_coord(anchor.primary_orf.transcript_start - 1)
            anchor_stop_codon = anchor.get_genome_coord_from_transcript_coord(anchor.primary_orf.transcript_stop - 1)
            
            if gencode:
                alt_transcripts = [tx for tx in chain(gc_transcripts.values(), pb_transcripts.values()) if tx is not anchor]
            else:
                alt_transcripts = [tx for tx in pb_transcripts.values() if tx is not anchor]
            
            for alternative in alt_transcripts:
                pblocks = ()
                with open(log_file, 'a') as flog:
                    with ExceptionLogger(info=f'{anchor}, {alternative}', output=flog):
                        other_start_codon = alternative.get_genome_coord_from_transcript_coord(alternative.primary_orf.transcript_start - 1)
                        other_stop_codon = alternative.get_genome_coord_from_transcript_coord(alternative.primary_orf.transcript_stop - 1)
                        anchor_starts_upstream = anchor_start_codon <= other_start_codon
                        anchor_stops_upstream = anchor_stop_codon <= other_stop_codon

                        alternative_length = alternative.protein.length
                        tx_aln = TranscriptAlignment.from_transcripts(anchor, alternative)
                        cd_aln = CodonAlignment.from_proteins(anchor.protein, alternative.protein)
                        pr_aln = ProteinAlignment.from_proteins(anchor.protein, alternative.protein)
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
                            'anchor': anchor.name,
                            'other': alternative.name,
                            'pblock_category': pblock.category.name,
                            'pblock_anchor_start': pblock.anchor_range.start,
                            'pblock_anchor_stop': pblock.anchor_range.stop,
                            'pblock_other_start': pblock.other_range.start,
                            'pblock_other_stop': pblock.other_range.stop,
                            'cblock_category': cblock.category.name,
                            'cblock_anchor_start': cblock.anchor_range.start,
                            'cblock_anchor_stop': cblock.anchor_range.stop,
                            'cblock_other_start': cblock.other_range.start,
                            'cblock_other_stop': cblock.other_range.stop,
                            'tblock_category': tblock.category.name if tblock else '',
                            'tblock_anchor_start': tblock.anchor_range.start if tblock else '',
                            'tblock_anchor_stop': tblock.anchor_range.stop if tblock else '',
                            'tblock_other_start': tblock.other_range.start if tblock else '',
                            'tblock_other_stop': tblock.other_range.stop if tblock else '',
                            'events': get_event_code(events),
                            'compound_splicing': any(set(events).intersection(compound_event.members) for compound_event in tx_aln.events if isinstance(compound_event, SpliceEvent) and len(compound_event.members) > 1),
                            'affects_up_start': anchor_starts_upstream and cblock is anchor_start_cblock or not anchor_starts_upstream and cblock is other_start_cblock,
                            'affects_down_start': anchor_starts_upstream and cblock is other_start_cblock or not anchor_starts_upstream and cblock is anchor_start_cblock,
                            'affects_up_stop': anchor_stops_upstream and cblock is anchor_stop_cblock or not anchor_stops_upstream and cblock is other_stop_cblock,
                            'affects_down_stop': anchor_stops_upstream and cblock is other_stop_cblock or not anchor_stops_upstream and cblock is anchor_stop_cblock,
                            'up_start_events': '',
                            'down_start_events': '',
                            'up_stop_events': '',
                            'down_stop_events': '',
                            # 'anchor_length': principal_length,
                            # 'other_length': alternative_length,
                            'cblock_anchor_seq': anchor.protein.sequence[cblock.anchor_range.start:cblock.anchor_range.stop],
                            'cblock_other_seq': alternative.protein.sequence[cblock.other_range.start:cblock.other_range.stop],                        
                        }
                        if cblock is anchor_start_cblock:
                            start_events = get_event_code(i.data for i in tx_aln.anchor_events.overlap(anchor.primary_orf.transcript_start - 1, anchor.primary_orf.transcript_start + 2) if isinstance(i.data, BasicTranscriptEvent))
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
                            stop_events = get_event_code(i.data for i in tx_aln.anchor_events.overlap(anchor.primary_orf.transcript_stop - 3, anchor.primary_orf.transcript_stop) if isinstance(i.data, BasicTranscriptEvent))
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

    print(f'Loading gene and transcript names for {chr}...')
    with db.get_session() as session:
        if gene_to_anchor_tx:
            condition = and_((Gene.chromosome_id == chr), Gene.name.in_(gene_to_anchor_tx))
        else:
            condition = (Gene.chromosome_id == chr)
        rows = session.execute(
                select(Gene.name, Transcript.name, Transcript.accession, Transcript.type).
                select_from(Protein).
                join(Protein.orf).
                join(ORF.transcript).
                join(Transcript.gene).
                where(condition)
        ).all()
        for gene_name, tx_name, tx_acc, tx_type in rows:
            gc_txs = gene_to_gc_transcripts.setdefault(gene_name, [])
            pb_txs = gene_to_pb_transcripts.setdefault(gene_name, [])
            if tx_type == 'gencodetranscript':
                gc_txs.append(tx_name)
            elif tx_type == 'pacbiotranscript':
                pb_txs.append(tx_acc)

    records = []
    # with mp.Pool() as p:
    t = tqdm(desc='Processing genes', total=len(gene_to_gc_transcripts), unit='gene', file=sys.stdout)
    for result in map(process_gene, gene_to_gc_transcripts.keys()):
        records.extend(result)
        t.update()
    chr_df = pd.DataFrame.from_records(records)
    return chr_df


def get_cblocks(db_name: str, output_dir: 'Path', log_dir: 'Path', gencode: bool, gene_to_anchor_tx: dict[str, str]):
    chrs = [f'chr{i}' for i in list(range(22, 23))]  # TODO: 
    dfs: dict[str, pd.DataFrame] = dict()

    for chr in chrs:
        df_file = output_dir/f'cblocks-{chr}.tsv'
        try:
            dfs[chr] = pd.read_csv(df_file, sep='\t')
        except:
            log_file = log_dir/f'{chr}.txt'
            dfs[chr] = process_chr(chr, db_name, log_file, gencode, gene_to_anchor_tx)
            dfs[chr].to_csv(df_file, sep='\t', index=False)
    
    cblock_df = pd.concat(dfs.values(), keys=dfs.keys(), names=['chr', 'row']).fillna(value='').reset_index().drop(columns='row')
    # cblock_df['other_accession'] = cblock_df['other'].str.split('|').str.get(1)
    return cblock_df


def get_pblocks(cblock_df: pd.DataFrame, output_dir: 'Path'):
    pblock_attrs = ['chr', 'anchor', 'other', 'pblock_category', 'pblock_anchor_start', 'pblock_anchor_stop']
    pblock_groups = cblock_df.groupby(pblock_attrs)

    pblocks = pd.DataFrame(index=pd.Index(data=pblock_groups.indices.keys(), name=pblock_attrs))

    # pblock_attrs = pblocks.index.to_frame()['pblock'].str.extract(r'\w\((\d+):(\d+)\|(\d+):(\d+)\)').astype(int)

    pblock_cat = pd.CategoricalDtype(['I', 'D', 'S'], ordered=True)

    pblocks['category'] = pblocks.index.to_frame()['pblock'].str.get(0).astype(pblock_cat)

    pblocks['length_change'] = pblocks['']
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

    pblocks.to_csv('', sep='\t')
