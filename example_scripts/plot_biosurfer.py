# %%

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from biosurfer.core.constants import APPRIS, CTerminalChange, NTerminalChange
from IPython.display import display
from matplotlib.patches import Patch, PathPatch
import os
import pickle


#%%
def run_plot():
    pblocks, output_dir = run_pickle()
    plot_hal(pblocks, output_dir)
#%%
def plot_hal(pblocks, output_dir):

    nterm_pblocks = pblocks[~pblocks['nterm'].isna() & (pblocks['nterm'] != NTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.ALTERNATIVE_ORF)].copy()
    nterm_pblocks['nterm'] = nterm_pblocks['nterm'].cat.remove_unused_categories()
    nterm_pblocks['altTSS'] = nterm_pblocks['events'].apply(lambda x: x.intersection('BbPp')).astype(bool)
    display(pd.crosstab(nterm_pblocks['up_start_cblock'], nterm_pblocks['down_start_cblock'], margins=True))

    nterm_palette = dict(zip(NTerminalChange, sns.color_palette('viridis_r', n_colors=5)[:-1]))

    nterm_fig = plt.figure(figsize=(3, 4))
    ax = sns.countplot(
        data = nterm_pblocks,
        y = 'nterm',
        order = (NTerminalChange.MUTUALLY_EXCLUSIVE, NTerminalChange.DOWNSTREAM_SHARED, NTerminalChange.UPSTREAM_SHARED, NTerminalChange.MUTUALLY_SHARED),
        palette = nterm_palette,
        linewidth = 0,
        saturation = 1,
    )
    ax.set(xlabel='Number of alternative isoforms', ylabel=None, yticklabels=[])
    plt.savefig(output_dir + '/nterm-class-counts.svg', dpi=200, facecolor=None)

    # %%
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
    plt.savefig(output_dir + '/nterm-length-change-dist.svg', dpi=200, facecolor=None)

    # %%
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
    plt.savefig(output_dir + '/nterm-altTSS-counts.svg', dpi=200, facecolor=None)

    # %%
    cterm_pblocks = pblocks[~pblocks['cterm'].isna() & (pblocks['nterm'] != NTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.ALTERNATIVE_ORF) & (pblocks['cterm'] != CTerminalChange.UNKNOWN)].copy()
    cterm_pblocks['cterm'] = cterm_pblocks['cterm'].cat.remove_unused_categories()
    cterm_pblocks['APA'] = cterm_pblocks['events'].apply(lambda x: x.intersection('BbPp')).astype(bool)

    display(pd.crosstab(cterm_pblocks['up_stop_cblock'], cterm_pblocks['down_stop_cblock'], margins=True))

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
    plt.savefig(output_dir + '/cterm-class-counts.svg', dpi=200, facecolor=None)

    # %%
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
    # cterm_splice_subcats = pd.DataFrame(
    #     {
    #         'EXIT (APA)': cterm_pblocks['up_stop_events'] == 'P',
    #         'EXIT (intron)': cterm_pblocks['up_stop_events'] == 'I',
    #         'EXIT (donor)': cterm_pblocks['up_stop_events'] == 'D',
    #         'ATE (multiple)': cterm_pblock_events.isin({('B', 'b'), ('b', 'B')}) & ~single_ATE,
    #         'ATE (single)': single_ATE,
    #         'poison exon': cterm_pblocks['up_stop_events'] == 'E',
    #         'other': [True for _ in cterm_pblocks.index]
    #     }
    # )
    cterm_pblocks['splice_subcat'] = cterm_splice_subcats.idxmax(axis=1).astype(pd.CategoricalDtype(cterm_splice_subcats.columns, ordered=True))

    cterm_splice_palette_dict = dict(zip(
        cterm_splice_subcats.columns,
        cterm_splice_palette[0:1] + cterm_splice_palette[1:2] + cterm_splice_palette[2:3] + ['#bbbbbb']
    ))

    # splice_subcat_order = cterm_pblocks[cterm_pblocks['cterm'] == CTerminalChange.SPLICING]['splice_subcat'].value_counts().index
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

    plt.savefig(output_dir + '/cterm-splicing-subcats.svg', dpi=200, facecolor=None, bbox_inches='tight')

    # %%
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

    plt.savefig(output_dir + '/cterm-frameshift-subcats.svg', dpi=200, facecolor=None, bbox_inches='tight')

    # %%
    cterm_event_counts = cterm_pblocks.groupby('cterm').events.value_counts().rename('count')
    cterm_examples = cterm_pblocks.groupby(['cterm', 'events']).sample(1, random_state=329).join(cterm_event_counts, on=['cterm', 'events']).sort_values(['cterm', 'count'], ascending=False).set_index(['cterm', 'events'])

    # %%
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
    plt.savefig(output_dir + '/internal-pblock-counts.svg', dpi=200, facecolor=None)

    # %%
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
    plt.savefig(output_dir + '/internal-pblock-events.svg', dpi=200, facecolor=None, bbox_inches='tight')

    # %%
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
    plt.savefig(output_dir + '/internal-pblock-events-split.svg', dpi=200, facecolor=None, bbox_inches='tight')

    # %%
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
    plt.savefig(output_dir + '/internal-pblock-compound-events.svg', dpi=200, facecolor=None, bbox_inches='tight')

def run_pickle():
    pickle_in = open("biosurfer.pickle","rb")
    pblocks = pickle.load(pickle_in)
    output_dir = pickle.load(pickle_in)
    return pblocks, output_dir


if __name__ == "__main__":
    run_plot()