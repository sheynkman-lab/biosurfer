import multiprocessing as mp
from operator import attrgetter

import matplotlib.pyplot as plt
from biosurfer.core.database import Database
from biosurfer.core.helpers import ExceptionLogger
from biosurfer.core.models.biomolecules import Gene
from biosurfer.plots.plotting import IsoformPlot
from tqdm import tqdm

db = Database('gencode')

#region
all_genes = (
    'A2ML1',
    'ACO2',
    'ACVRL1',
    'ADGRG1',
    'ANKUB1',
    'APOOL',
    'ARL13B',
    'ASPG',
    'BCL10',
    'BNIP2',
    'CALCOCO2',
    'CAPN3',
    'CCDC88A',
    'CCM2',
    'CMTM5',
    'CYP4F11',
    'DENND1B',
    'DIAPH1',
    'DMRTC2',
    'DNAJC24',
    'EIF6',
    'ENTPD8',
    'EVPL',
    'FBXL16',
    'FKBPL',
    'FYB1',
    'GAS2L1',
    'GATD3A',
    'GIMAP5',
    'GNG5P2',
    'GSDMC',
    'GTF2H3',
    'IL31',
    'ITGAD',
    'JOSD2',
    'KIAA0825',
    'KRTAP10-3',
    'LY6D',
    'LYPD1',
    'LZTS1',
    'MAST4',
    'MOCOS',
    'MOGAT1',
    'MORC2',
    'MUCL3',
    'MYO18A',
    'NANOS3',
    'NAP1L5',
    'NOLC1',
    'NPEPPS',
    'NR0B2',
    'NTPCR',
    'NXPH2',
    'OR10H4',
    'OR2A4',
    'OR52E8',
    'OR52N1',
    'PARVB',
    'PDAP1',
    'PDE6G',
    'PFDN6',
    'PLEKHA7',
    'PON3',
    'PPAT',
    'PRH1-PRR4',
    'PRR9',
    'PRSS51',
    'PRYP3',
    'PTPRF',
    'RANBP3L',
    'RHBDD1',
    'RHOD',
    'RSBN1L',
    'SCN11A',
    'SEMA4A',
    'SGK2',
    'SLC24A1',
    'SLC35G2',
    'SLX4',
    'SPRR2F',
    'SRCAP',
    'STAT5B',
    'SULT2B1',
    'TBC1D8',
    'TEAD1',
    'TEX101',
    'TGIF2',
    'TLE3',
    'TLX2',
    'TP53INP2',
    'TP53TG5',
    'TPGS2',
    'UBE2U',
    'UFSP2',
    'UMPS',
    'UROS',
    'WDR5B',
    'WNT8A',
    'ZNF487',
    'ZSWIM2',
)
#endregion

def process(gene_name):
    with db.get_session() as session:
        gene = Gene.from_name(session, gene_name)
        with ExceptionLogger(gene_name):
            isoplot = IsoformPlot(sorted(gene.transcripts, key=attrgetter('appris'), reverse=True))
            isoplot.draw_all_isoforms()
            isoplot.draw_features()
            isoplot.draw_frameshifts()
            isoplot.draw_legend()
            isoplot.fig.set_size_inches(12, 0.5 + 0.5*len(gene.transcripts))
            plt.savefig(f'../output/test/{gene_name}.png', dpi=150, facecolor='w', bbox_inches='tight')
            plt.close(isoplot.fig)

input = all_genes

with mp.Pool() as p:
    _ = list(tqdm(p.imap(process, input), desc='Processing genes', total=len(input)))
