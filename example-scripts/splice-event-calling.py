# %%
import os
from time import time

import matplotlib.pyplot as plt
import pandas as pd
from biosurfer.core.alignments_old import Alignment as ProteinAlignmentOld
from biosurfer.core.alignments import CodonAlignment, TranscriptAlignment, ProteinAlignment
from biosurfer.core.constants import \
    CodonAlignmentCategory as CodonAlignCat
from biosurfer.core.database import Database
from biosurfer.core.models.biomolecules import Transcript
from biosurfer.core.splice_events import (AcceptorSpliceEvent,
                                          DonorSpliceEvent, ExonBypassEvent,
                                          ExonSpliceEvent, IntronSpliceEvent)
from biosurfer.plots.plotting import IsoformPlot
from tqdm import tqdm

df = pd.read_csv('complex-splice-events.csv', names=['anchor', 'other', 'score'])

db = Database('gencode')
session = db.get_session()
txs = Transcript.from_names(session, pd.concat([df['anchor'], df['other']]).unique())

# %%
all_alns: dict[tuple[str, str], tuple['TranscriptAlignment', 'CodonAlignment', 'ProteinAlignmentOld']] = dict()
t0 = time()
for a, o, _ in df.sort_values('anchor').itertuples(index=False):
    anchor: 'Transcript' = txs[a]
    other: 'Transcript' = txs[o]

    try:
        pr_aln = ProteinAlignment.from_proteins(anchor.protein, other.protein)
    except ValueError as e:
        pr_aln = None
        print(f'{a}, {o}: {e}')
    cd_aln = CodonAlignment.from_proteins(anchor.protein, other.protein)
    tx_aln = TranscriptAlignment.from_transcripts(anchor, other)
    pr_aln_old = ProteinAlignmentOld(anchor.protein, other.protein)
    all_alns[a, o] = tx_aln, cd_aln, pr_aln, pr_aln_old
t1 = time()

# %%
for (a, o), (tx_aln, cd_aln, pr_aln, pr_aln_old) in all_alns.items():
    print(tx_aln)
    for event in tx_aln.events:
        print(f'\t{type(event).__name__}')
        for e in event.members:
            print(f'\t\t{e}')
    print(f'{a} events')
    for i in sorted(tx_aln.anchor_events):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t' + getattr(i.data, 'code', type(i.data).__name__))
    print(f'{o} events')
    for i in sorted(tx_aln.other_events):
        print(f'\t{i.begin:5d}\t{i.end:5d}\t' + getattr(i.data, 'code', type(i.data).__name__))
    print(f'{cd_aln} blocks')
    for block in cd_aln.blocks:
        print(f'\t{block.category}\t{block.anchor_range.start:5d} {block.anchor_range.stop:5d}\t| {block.other_range.start:5d} {block.other_range.stop:5d}')
    print(f'{pr_aln_old} cblocks (old)')
    for block in pr_aln_old.transcript_blocks:
        anchor_range = (block.anchor_residues[0].position - 1, block.anchor_residues[-1].position) if block.anchor_residues else (-1, -1)
        other_range = (block.other_residues[0].position - 1, block.other_residues[-1].position) if block.other_residues else (-1, -1)
        print(f'\t{block.category}\t{anchor_range[0]:5d} {anchor_range[1]:5d}\t| {other_range[0]:5d} {other_range[1]:5d}')
    print(f'{pr_aln} blocks')
    for block in pr_aln.blocks:
        print(f'\t{block.category}\t{block.anchor_range.start:5d} {block.anchor_range.stop:5d}\t| {block.other_range.start:5d} {block.other_range.stop:5d}')
    print(f'{pr_aln_old} pblocks (old)')
    for block in pr_aln_old.protein_blocks:
        anchor_range = (block.anchor_residues[0].position - 1, block.anchor_residues[-1].position) if block.anchor_residues else (-1, -1)
        other_range = (block.other_residues[0].position - 1, block.other_residues[-1].position) if block.other_residues else (-1, -1)
        print(f'\t{block.category}\t{anchor_range[0]:5d} {anchor_range[1]:5d}\t| {other_range[0]:5d} {other_range[1]:5d}')

# %%
EVENT_COLORS = {
    IntronSpliceEvent: '#e69138',
    DonorSpliceEvent: '#6aa84f',
    AcceptorSpliceEvent: '#674ea7',
    ExonSpliceEvent: '#3d85c6',
    ExonBypassEvent: '#bebebe',
}

BLOCK_COLORS = {
    CodonAlignCat.TRANSLATED: '#9bf3ff',
    CodonAlignCat.INSERTION: '#05e0ff',
    CodonAlignCat.FRAME_AHEAD: '#fff099',
    CodonAlignCat.FRAME_BEHIND: '#ffd700',
    CodonAlignCat.UNTRANSLATED: '#ff99ce',
    CodonAlignCat.DELETION: '#ff0082',
    CodonAlignCat.EDGE: '#8270c1',
    CodonAlignCat.COMPLEX: '#aaaaaa'
}

force_plotting = False
for a, others in df.groupby('anchor')['other']:
    anchor = txs[a]
    transcripts = [anchor] + [txs[o] for o in others]
    fig_path1 = f'../output/alignment-redesign/splice-events/{txs[a].gene.name}.png'
    fig_path2 = f'../output/alignment-redesign/cblocks/{txs[a].gene.name}.png'
    fig_path1_missing = not os.path.isfile(fig_path1)
    fig_path2_missing = not os.path.isfile(fig_path2)
    if force_plotting or fig_path1_missing:
        isoplot = IsoformPlot(transcripts)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()

        for i, other in enumerate(others, start=1):
            for event in all_alns[a, other][0].events:
                for base_event in event.members:
                    isoplot.draw_region(
                        i,
                        base_event.start.coordinate,
                        base_event.stop.coordinate,
                        y_offset = -0.5,
                        height = 0.1,
                        facecolor = EVENT_COLORS[type(base_event)]
                    )
        isoplot.fig.set_size_inches(10, 0.2 + 0.4 * len(transcripts))
        plt.savefig(fig_path1, dpi=200, bbox_inches='tight', facecolor='w')
        print('saved '+fig_path1)
        plt.close(isoplot.fig)
    if force_plotting or fig_path2_missing:
        isoplot = IsoformPlot(transcripts)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        for i, o in enumerate(others, start=1):
            isoplot.draw_codon_alignment_blocks(all_alns[a, o][1])
        isoplot.fig.set_size_inches(10, 0.2 + 0.4 * len(transcripts))
        plt.savefig(fig_path2, dpi=200, bbox_inches='tight', facecolor='w')
        print('saved '+fig_path2)
        plt.close(isoplot.fig)

# %%
print(f'{t1 - t0:.3g} s total')
print(f'{(t1 - t0)/len(all_alns):.3g}s per alignment')

# %%
