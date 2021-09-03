from operator import attrgetter
from os.path import isfile
from warnings import filterwarnings

import matplotlib.pyplot as plt
import pandas as pd
from biosurfer.analysis.sqtl import get_pblocks_related_to_junction
from biosurfer.core.alignments import export_annotated_pblocks_to_tsv
from biosurfer.core.constants import ProteinLevelAlignmentCategory
from biosurfer.core.models import Chromosome, Gene, Junction
from biosurfer.plots.plotting import IsoformPlot
from IPython.display import display
from matplotlib._api.deprecation import MatplotlibDeprecationWarning

filterwarnings("ignore", category=MatplotlibDeprecationWarning)

colors = {
    ProteinLevelAlignmentCategory.DELETION: 'mediumvioletred',
    ProteinLevelAlignmentCategory.INSERTION: 'goldenrod',
    ProteinLevelAlignmentCategory.SUBSTITUTION: 'lightseagreen'
}

data_dir = '../data/bone'
output_dir = '../output/bone'
sqtls = pd.read_csv(f'{data_dir}/fibroblast_coloc_sqtls.csv')
sqtls = sqtls[sqtls['gene_type'] == 'protein_coding']
chromosomes = Chromosome.from_names(set(sqtls['chr']))
genes = Gene.from_names(set(sqtls['gene_name']))

for index, row in sqtls.iterrows():
    gene = genes[row['gene_name']]
    junc = Junction(row['event_start'], row['event_end'], chromosomes[row['chr']], gene.strand)
    
    pblocks, using, not_using = get_pblocks_related_to_junction(junc, gene.transcripts)
    if not all((using, not_using)):
        continue
    print(junc)
    print(f'\ttranscripts using: {using}')
    print(f'\ttranscripts not using: {not_using}')
    export_annotated_pblocks_to_tsv(f'{output_dir}/{gene.name}_{junc.donor}_{junc.acceptor}.tsv', pblocks)

    fig_path = f'{output_dir}/{gene.name}_{junc.donor}_{junc.acceptor}.png'
    if not isfile(fig_path):
        isoplot = IsoformPlot(sorted(using, key=attrgetter('appris')) + sorted(not_using, key=attrgetter('appris')))
        isoplot.draw_all_isoforms()
        isoplot.draw_frameshifts()
        isoplot.draw_background_rect(start=junc.donor, stop=junc.acceptor, facecolor='#ffffb7')

        for pblock in pblocks:
            anchor = isoplot.transcripts.index(pblock.parent.anchor.transcript)
            other = isoplot.transcripts.index(pblock.parent.other.transcript)
            if pblock.category is ProteinLevelAlignmentCategory.DELETION:
                start = pblock.anchor_residues[0].codon[1].coordinate
                stop = pblock.anchor_residues[-1].codon[1].coordinate
            else:
                start = pblock.other_residues[0].codon[1].coordinate
                stop = pblock.other_residues[-1].codon[1].coordinate
            isoplot.draw_region(
                other,
                start = start,
                stop = stop,
                y_offset = 0.5*isoplot.opts.max_track_width,
                height = 0.4*isoplot.opts.max_track_width,
                edgecolor = 'none',
                facecolor = colors[pblock.category]
            )

        isoplot.fig.set_size_inches(9, 0.5*len(gene.transcripts))
        plt.savefig(fig_path, facecolor='w', transparent=False, dpi=300, bbox_inches='tight')
        print('\tsaved '+fig_path)
        plt.close(isoplot.fig)
