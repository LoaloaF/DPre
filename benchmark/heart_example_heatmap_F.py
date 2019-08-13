import pandas as pd
import numpy as np
import os
import DPre

from mixed_sample_ident_benchmark import benchmark_multi

import matplotlib.pyplot as plt
import matplotlib as mpl
# os.chdir('benchmark')

# the target map assigns the closest matching target to each sample 
# (ideally, this target is rated as the most transcriptionally similar) 
trg_map = pd.read_csv('mappings/target_mapping.tsv', sep='\t', 
                      index_col=0, squeeze=True).dropna()
h_expr = ['expr/heart_expression.tsv', 
          'expr/hsliver_expression.tsv', 
          'expr/microglia_macro_corticon_expression.tsv']
# merge human expr to one expression to pool from when making synthetic combinations
h_expr = pd.concat([pd.read_csv(expr, sep='\t', index_col='ensg') 
                    for expr in h_expr], axis=1, sort=False)
heart_example = pd.DataFrame(['mean_hESC', 'mean_Fetal Heart', 'invitro DA neurons', np.nan], 
                             index=['start', 'target', 'sample2', 'sample3']).T

# define the proportions in which combinations are formed
bins = 19
# props of the target of interest
trg_props = np.arange(.05, 1, .05)
# props of the matching base (control/ null-state, of ESCs) for combs with 3 and 4 cell types (2 and 3 somatic noise)
comb_base23_props = np.append(np.arange(.9, .2, -.1), np.arange(.3, 0, -.025))
# props of the ct representing somatic noise in a 3-ct combination
comb2_trg2_props = np.append(np.arange(.05, .4, .05), np.arange(.3, 0, -.025))
# props of the 2 cts representing somatic noise in a 4-ct combination
comb3_trg23_props = np.append(np.arange(.025, .2, .025), np.arange(.15, 0, -.0125))
# glue together
props_arr = np.stack([comb_base23_props, comb2_trg2_props, trg_props], axis=1)
props_arr = np.vstack( [[1,0,0], [0,1,0], props_arr, [0, 0, 1]])
print(props_arr)
print(' Sum of rows:\n', props_arr.sum(axis=1)) # should be 1 for all

def do_heatmap(species, metric):
    res = benchmark_multi(None, h_expr, m_combs=None, h_combs=heart_example, 
                          metric=metric, return_single=True)
    t, smp, trg_cts = res
    trg_cts = np.concatenate([[trg_cts[0]], trg_cts[2:], [trg_cts[1]]])
    lbls = [name[:-5] for name in smp.names[:len(trg_cts)]]
    lbls = ['hESC', '$\mathit{in}$ $\mathit{vitro}$' + ' dopaminergic neurons', 'Fetal Heart']

    props_arr = np.stack([comb_base23_props, comb2_trg2_props, trg_props], axis=1)
    props_arr = np.vstack([[1,0,0], [0,1,0], props_arr, [0, 0, 1]])
    
    smp.reorder(list(reversed(smp.names)))
    sim = t._get_similarity(smp, metric)[0].xs('mean',1,0)

    scores = []
    res = []
    def rank_scores(row):
        row = row.sort_values(ascending=False)
        for trg_name in trg_cts:
            scores.append(row.index.get_loc(trg_name)+1)
        res.append(pd.Series(scores, index=trg_cts, name=row.name))
        scores.clear()
    sim.apply(rank_scores, axis=1)
    ranks = pd.concat(res, axis=1).T.reindex(trg_cts, axis=1)

    t = t.slice_elements(trg_cts)
    hm_width = 1.5
    pivot = True
    hm = t.target_similarity_heatmap(smp, 
                                metric=metric, 
                                differential=False,
                                plt_show=False, 
                                pivot=pivot,
                                heatmap_width = hm_width,
                                heatmap_height=.45,
                                heatmap_range=[270, 0],
                                hide_samplelabels = True,
                                show_samples_colorbar = True,
                                hide_distance_bar = True,
                                # samplelabels_space= .5,
                                HM_RIGHT = .5,
                                # HM_WSPACE = .5,
                                HM_Y_COLORBAR = .11 *len(trg_cts) *hm_width,
                                title = 'Accuracy of identification\non Fetal heart in mixed samples',
                                filename=None)
    fig = hm[0]
    axes = hm[1]
    sim = hm[2][0]
    custom_prop_hm = [ranks, props_arr, lbls]
    width, height = fig.get_size_inches()
    xlbl = trg_cts
    # ADJUST NORMAL HEATMAP
    sim = custom_prop_hm[0]
    perc_thresh = (sim.iloc(1)[-1]<(301*.05)).sum()
    req_prop = float(sim.iloc(1)[-1][sim.iloc(1)[-1]<(301*.05)].index[-1][-4:])

    ax = axes[2, 1]
    ax.hlines(perc_thresh-.5, -3.8, 3, linewidth=.6, clip_on=False)
    ax.annotate('{}%'.format(int(req_prop*100)), xy=(3.2, perc_thresh-1.7), annotation_clip=False,
                            rotation_mode='anchor', ha='center', va='center', zorder=200)

    [ax.text(2, y+.25, '*', ha="center", va="center", color="w", fontsize=4) 
        for y in np.arange(perc_thresh-1, -.1, -1)]
    
    hm_args = {'cmap': 'RdBu', 'vmin': 0, 'vmax': 270}
    im = ax.imshow(sim.values, aspect='auto', **hm_args)
    
    ax.set_xticks(np.arange(sim.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(sim.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="k", linestyle='-', linewidth=.5, alpha=.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_yticks(np.arange(0, sim.shape[0]))

    # PROPORTION HEATMAP
    import matplotlib as mpl
    from matplotlib.colors import LinearSegmentedColormap
    cm = LinearSegmentedColormap.from_list('cm', 
                                        ['w', 'k'], 
                                        9)
    ax = axes[2,0]
    prop_data = custom_prop_hm[1]
    # prop_data = prop_data[::-1,:]
    [spine.set_visible(True) for spine in ax.spines.values()]
    # The first column in this heatmap is somehow broken. DOesn't want to draw in left column???
    print(prop_data)
    im = ax.imshow(prop_data, cmap=cm, aspect='auto')
    ax.tick_params(labelleft=False, labelbottom=True) 
    ax.set_ylim(0, 20)

    fs = DPre.config.FONTS
    axes[3,1].set_xticklabels(xlbl, rotation=30, ha='right', va='top', 
                        fontsize=fs, rotation_mode='anchor', x=.5)
    ax.set_xticks([0,1,2])
    ax.set_xticklabels(custom_prop_hm[2], rotation=30, ha='right', va='top', 
                        fontsize=fs, rotation_mode='anchor', y=-.053)
        
    ax.set_xticks(np.arange(prop_data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(prop_data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="k", linestyle='-', linewidth=.5, alpha=.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    ax.set_yticks(np.arange(0, prop_data.shape[0]))
    ax.tick_params(left=True)

    at = ((DPre.config.CB_LEFT_SEC)/width, 1- (DPre.config.CB_TOP)/height, 
        DPre.config.CB_WIDTH/width, DPre.config.CB_HEIGHT/height)
    cax = fig.add_axes(at)
    cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
    
    xts = np.arange(.1, .91, .4)
    cb.set_ticks(xts)
    cb.ax.set_xticklabels([str(int(xt*100))+'%' for xt in xts])
    if pivot:
        cb.ax.tick_params(labelrotation=90)
    cb.ax.set_xlabel('Proportion in\nsynthetic sample')
    cb.ax.get_xaxis().set_label_position('top')
    fig.savefig('plots/heart_example_heatmap.svg')
    fig.savefig('plots/heart_example_heatmap.png')
do_heatmap(species='human', metric='cosine')