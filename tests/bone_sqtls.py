import os
from itertools import groupby, product
from operator import attrgetter
from statistics import median
from warnings import filterwarnings

import matplotlib.pyplot as plt
import pandas as pd
from biosurfer.analysis.sqtl import (get_event_counts,
                                     get_pblocks_related_to_junction,
                                     junction_causes_knockdown_in_pair,
                                     split_transcripts_on_junction_usage)
from biosurfer.core.alignments import (TranscriptBasedAlignment,
                                       export_annotated_pblocks_to_tsv)
from biosurfer.core.models import Chromosome, Gene, Junction
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

data_dir = '../data/bone'
output_dir = '../output/bone'
sqtls = pd.read_csv(f'{data_dir}/fibroblast_coloc_sqtls.csv')
sqtls = sqtls[sqtls['gene_type'] == 'protein_coding']
chromosomes = Chromosome.from_names(set(sqtls['chr']))
genes = Gene.from_names(set(sqtls['gene_name']))

records = []
for index, row in sqtls.iterrows():
    gene = genes[row['gene_name']]
    junc = Junction(row['event_start'], row['event_end'], chromosomes[row['chr']], gene.strand)
    junc_info = {
        'chr': row['chr'],
        'strand': str(junc.strand),
        'donor': row['event_start'],
        'acceptor': row['event_end'],
        'gene': gene.name,
        'H4': row['H4'],
        'pairs': 0
    }

    using, not_using = split_transcripts_on_junction_usage(junc, gene.transcripts)
    using = sorted(using, key=attrgetter('appris'))
    not_using = sorted(not_using, key=attrgetter('appris'))
    if not all((using, not_using)):
        continue
    print(junc)
    print(f'\t{len(using)} transcripts using: {using}')
    print(f'\t{len(not_using)} transcripts not using: {not_using}')
    
    pairs = list(product(not_using, using))
    n_pairs = len(pairs)
    alns = [TranscriptBasedAlignment(tx1.protein, tx2.protein) for tx1, tx2 in pairs]
    junc_info['pairs'] = n_pairs

    pblocks = get_pblocks_related_to_junction(junc, alns)

    # calculate fraction of pairs where junction-lacking isoform is not NMD but junction-using isoform is NMD
    junc_info['NMD'] = sum(tx2.primary_orf.nmd and not tx1.primary_orf.nmd for tx1, tx2 in pairs) / n_pairs

    # calculate median change in sequence length for junction-related pblocks
    junc_info['median_delta_length'] = median(pblock.delta_length for pblock in pblocks) if pblocks else 0.0

    # classify knockdown vs. alteration for each pair
    # threshold = -not_using[0].protein.length * 2 // 5
    pair_pblocks = groupby(pblocks, key=attrgetter('anchor', 'other'))
    knockdown_pair_count = sum(junction_causes_knockdown_in_pair(pblocks=list(pblock_group)) for _, pblock_group in pair_pblocks)
    junc_info['knockdown_frequency'] = knockdown_pair_count / n_pairs

    counts = get_event_counts(pblocks)
    freqs = {event.name: count/n_pairs for event, count in counts.items()}
    junc_info.update(freqs)

    records.append(junc_info)

    if not os.path.isdir(f'{output_dir}/{gene.name}'):
        os.mkdir(f'{output_dir}/{gene.name}')
    export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks)

    fig_path = f'{output_dir}/{gene.name}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
    if not os.path.isfile(fig_path):
        isoplot = IsoformPlot(using + not_using)
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#ffffb7')
        for pblock in pblocks:
            isoplot.draw_protein_block(pblock)
        isoplot.fig.set_size_inches(9, 0.5*len(gene.transcripts))
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
        print('\tsaved '+fig_path)
        plt.close(isoplot.fig)

sqtls_augmented = pd.DataFrame.from_records(records)
display(sqtls_augmented)
sqtls_augmented.to_csv(f'{output_dir}/fibroblast_sqtls_annotated.tsv', sep='\t')
