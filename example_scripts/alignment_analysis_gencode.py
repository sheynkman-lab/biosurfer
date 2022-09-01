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
        #print('\n Database ' + db_name + ' does not exist. Pleaase create a new database or use an existing database ... \n')
        raise ValueError('\n Database ' + db_name + ' does not exist. Please create a new database or use an existing one ... \n')

   
def run_hybrid_alignment(db_name, output_path):
    """ Runs hybrid alignment on input database.
    Args:
        db_name: User input database name
        output_path: Directory to write output to.

    Returns:
      Nothing
    """

    # Sanity check for the directory name and database name
    db_check(db_name)


    plt.rcParams['svg.fonttype'] = 'none'
    sns.set_style('whitegrid')

    db_name = db_name
    db = Database(db_name)

    output_dir = Path(output_path)
    output_dir = output_dir.resolve()

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
                            'pblock': str(pblock),
                            #'pblock_anchor_start': pblock.anchor_range.start,
                            #'pblock_anchor_stop': pblock.anchor_range.stop,
                            'cblock': str(cblock),
                            'cblock category': cblock.category,
                            'cblock_anchor_start': cblock.anchor_range.start,
                            'cblock_anchor_stop': cblock.anchor_range.stop,
                            'cblock_other_start': cblock.other_range.start,
                            'cblock_other_stop': cblock.other_range.stop,
                            'tblock': str(tblock),
                            #'tblock_anchor_start': tblock.anchor_range.start,
                            #'tblock_anchor_stop': tblock.anchor_range.stop,
                            'events': get_event_code(events),
                            'compound_splicing': any(set(events).intersection(compound_event.members) for compound_event in tx_aln.events if isinstance(compound_event, SpliceEvent) and len(compound_event.members) > 1),
                            'anchor_aa': principal.protein.sequence[cblock.anchor_range.start:cblock.anchor_range.stop],
                            #'tblock_anchor_seq': principal.protein.sequence[tblock.anchor_range.start:tblock.anchor_range.stop],
                            #'pblock_anchor_seq': principal.protein.sequence[pblock.anchor_range.start:pblock.anchor_range.stop],
                            'other_aa': alternative.protein.sequence[cblock.other_range.start:cblock.other_range.stop],
                            #'tblock_anchor_seq': alternative.protein.sequence[tblock.other_range.start:tblock.other_range.stop],
                            #'pblock_anchor_seq': alternative.protein.sequence[pblock.other_range.start:pblock.other_range.stop],
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
            # fig_path = output_dir/f'chr19/{gene_name}.svg'
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

        df_path = output_dir/f'alignment-analysis-{chr}.tsv'
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
    display(df)

    # sys.exit()

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

    # %%
    # changed s type to string to print in text files.
    def get_genes(pblocks_slice: 'pd.DataFrame'):
        return pblocks_slice.reset_index()['anchor'].str.extract(r'(\w+)-\d+', expand=False).unique()
                 
    
    with open(output_dir/f'genes-mxs.txt', 'w') as f:
        f.writelines(str(s)+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.MUTUALLY_EXCLUSIVE]))

    with open(output_dir/f'genes-sds.txt', 'w') as f:
        f.writelines(str(s)+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.DOWNSTREAM_SHARED]))

    with open(output_dir/f'genes-sus.txt', 'w') as f:
        f.writelines(str(s)+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.UPSTREAM_SHARED]))

    with open(output_dir/f'genes-mss.txt', 'w') as f:
        f.writelines(str(s)+'\n' for s in get_genes(pblocks[pblocks['nterm'] == NTerminalChange.MUTUALLY_SHARED]))

    # %%
    display(pblocks)

    # Pickle the pblocks datafrom for vis.
    pickle_out = open("biosurfer.pickle","wb")
    pickle.dump(pblocks, pickle_out)
    pickle.dump(output_path, pickle_out)
    pickle_out.close()


if __name__ == '__main__':
    run_hybrid_alignment(db_name, output_path)

    
