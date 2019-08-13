import pandas as pd
import numpy as np
import os
import DPre

import matplotlib.pyplot as plt
import matplotlib as mpl
os.chdir('benchmark')

"""Boxplots"""
def make_boxplot(which, metric, fname, single_benchmark_data=None):
    # function that returns a collection of symetric x locs according to the 
    # length and base of the input x loc 
    def shift_x(data, postions, x_stepsize=.05):
        xs = []
        for d, p in zip(data, postions):
            x = np.array([])
            for v in np.unique(d, return_counts=True)[1]:
                if v == 1:
                    x = np.append(x, [p])
                elif v == 2:
                    x = np.append(x, [p-x_stepsize/2, p+x_stepsize/2])
                elif v %2 == 0:
                    arr = np.linspace(p -(v/2)*x_stepsize, p +(v/2)*x_stepsize, v)
                    x = np.append(x, arr)
                elif v %2 == 1:
                    arr = np.linspace(p -((v-1)/2)*x_stepsize, p +((v-1)/2)*x_stepsize, v)
                    x = np.append(x, arr)
            xs.append(x)
        return np.concatenate(xs)
        
    # SINGLE DATA
    if which == 'ct_ident':
        sing_d = pd.read_csv(fname, sep='\t')
        sing_d_euc = sing_d[sing_d['metric'] == 'euclid']['rank'].values
        sing_d_pea = sing_d[sing_d['metric'] == 'pearson']['rank'].values
        sing_d_cos = sing_d[sing_d['metric'] == 'cosine']['rank'].values
        sing_d_int = sing_d[sing_d['metric'] == 'intersect']['rank'].values
        sing_d_euc.sort()
        sing_d_pea.sort()
        sing_d_cos.sort()
        sing_d_int.sort()
        data = [sing_d_euc, sing_d_pea, sing_d_cos, sing_d_int]
        
        xts = np.arange(0, 1.8, .5)
        xts[-1] += .3
        fig = plt.figure(figsize=(2.2, 2.2))

    # MULTI DATA
    elif which == 'mixed_samples':
        valid_samples = pd.read_csv(single_benchmark_data, sep='\t', index_col=(0, 1)).xs(metric)
        mouse = valid_samples.iloc[7:].index[valid_samples['rank'].iloc[7:] <14]
        human = valid_samples.iloc[:7].index[valid_samples['rank'].iloc[:7] <16]
        valid_samples = human.append(mouse)

        multi_d = pd.read_csv(fname, sep='\t', index_col=0)
        all_samples = np.unique(multi_d['target'].values)
        inval = pd.Index(all_samples).drop(valid_samples, 'ignore').values

        # Samples that are not within 5% on their own. Drop.
        multi_d = multi_d[~(np.isin(multi_d['target'], inval) & (multi_d['required_prop'] >1) )]#& multi_d['sample2'].isna())]
        # Samples that require slightly more than 100% proportion due to merged datasets effects
        # plotted as 100$ required. Set to 1
        multi_d.required_prop.mask(multi_d.required_prop >1, 1.0, inplace=True)
        # sample combinations where 0% of the sample already score in 5%. Drop.
        multi_d = multi_d[multi_d.required_prop != 0]

        multi_d_n2 = multi_d[multi_d['sample2'].isna()]['required_prop'].values
        multi_d_n3 = multi_d[multi_d['sample2'].notna() & multi_d['sample3'].isna()]['required_prop'].values
        multi_d_n4 = multi_d[multi_d['sample2'].notna() & multi_d['sample3'].notna()]['required_prop'].values
        print('2Cts n: {}'.format(len(multi_d_n2)))
        print('3Cts n: {}'.format(len(multi_d_n3)))
        print('4Cts n: {}'.format(len(multi_d_n4)))
        
        multi_d_n2.sort()
        multi_d_n3.sort()
        multi_d_n4.sort()
        data = [multi_d_n2, multi_d_n3, multi_d_n4]
        xts = np.arange(0, 1.3, .5)
        fig = plt.figure(figsize=(1.6, 2.2))



    ax = fig.add_subplot(111)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    boxprops = dict(linewidth=.5)
    whiskerprops = dict(linewidth=.5)
    capprops = dict(linewidth=.5)
    flierprops = dict(linewidth=.1, markersize=.5, marker='^')
    medianprops = dict(linewidth=.5, color='k')
    width=.4
    pb = ax.boxplot(data, positions=xts,
                widths=width, showfliers=False,
                boxprops=boxprops, flierprops=flierprops, medianprops=medianprops, whiskerprops=whiskerprops, capprops=capprops)

    x = shift_x(data, postions=xts, x_stepsize=.04)
    y = np.concatenate(data)
    ax.scatter(x, y, s=1, zorder=200, color='#59a48c', edgecolor=None)

    # # SINGLE DATA
    if which == 'ct_ident':
        ax.set_xticklabels(['Euclidean distance', 'Pearson correlation', 'Cosine similarity', 'Intersect metric'], fontsize=6, 
                            rotation=30, ha='right', va='top')
        fig.subplots_adjust(bottom=.30, left=.32, right=.7)
        
        ax.set_xticks(xts)
        ax.set_title('Cell type identification\nbenchmark', fontsize=DPre.config.FONTS)
        ax.set_ylabel('Similarity rank', fontsize=DPre.config.FONTS, rotation=270, va='center', labelpad=6, ha='center')
        ax.spines['right'].set_visible(True)
        ax_r = ax.twinx()
        ax_r.tick_params(right=False, labelright=False)
        ax_r.set_ylabel('Similarity increase rank', fontsize=DPre.config.FONTS, rotation=270, labelpad=6, va='center', ha='center')
        ax_r.spines['top'].set_visible(False)
        ax.vlines(1.4, -10, 310, linewidth = 4, color='w', clip_on=False)

        ax.tick_params(bottom=False)
        yticks = np.arange(0,260, 50)
        yticks[0] = 1
        yticklbls = [str(t)+'.'for t in yticks]
        ax.set_yticks(yticks)
        ax.set_ylim((0, 300))
        ax.set_xlim(xts[0]-width/2 -.1, xts[-1]+width/2 +.1 +.3)
        ax.set_yticklabels(yticklbls, va='center', fontsize=DPre.config.FONTS)
        ax.yaxis.grid(linestyle='dashed', alpha=.6)
        ax.set_ylim((-1, ax.get_ylim()[1]))
        fig.savefig('plots/single_benchmark.png')
        fig.savefig('plots/single_benchmark.svg')
    
    # MULTI DATA
    elif which == 'mixed_samples':
        fig.subplots_adjust(bottom=.33, left=.38, right=.7)
        lbl = ('Minimum required proportion\nfor correct cell type\n'
               'identification in the top 5%')
        tit  = 'Accuracy of identification\nin mixed samples' 
        ax.set_title(tit, fontsize=DPre.config.FONTS)
        ax.set_ylabel(lbl, fontsize=DPre.config.FONTS-1)
        ax.set_xticklabels(['2 cell types', '3 cell types', '4 cell types'], 
                        fontsize=DPre.config.FONTS, rotation=30, ha='right')
        ax.yaxis.grid(linestyle='dashed', alpha=.6)
        yts = np.arange(0, 1.1, .25)
        ax.set_xlim(xts[0]-width/2 -.1, xts[-1]+width/2 +.1)
        ax.set_yticks(yts)
        ax.set_ylim(-.05, 1.05)
        ax.set_yticklabels([str(int(t*100)) +'%' for t in yts], fontsize=DPre.config.FONTS-1)
        
        # plt.show()
        fig.savefig('plots/{}_mixed_samples_benchmark.svg'.format(metric))
        fig.savefig('plots/{}_mixed_samples_benchmark.png'.format(metric))

metric = 'cosine'
cti_fname = 'benchmark_scores/ct_ident_benchmark_scores.tsv'
cti_which = 'ct_ident'
make_boxplot(cti_which, metric, cti_fname)

ms_fname = 'benchmark_scores/{}_mixed_sample_benchmark_scores_merged.tsv'.format(metric)
ms_which = 'mixed_samples'
make_boxplot(ms_which, metric, ms_fname, single_benchmark_data=cti_fname)