import pandas as pd
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import distance, linkage, dendrogram, optimal_leaf_ordering
from scipy.spatial.distance import pdist

import DPre.main.config as config
from DPre.main._logger import logger, spacer



def _add_mg_types(expr, down):
        orig_order = expr.columns.unique(0)
        updown_idx = ['up']*expr.shape[1]
        if down:
            updown_idx.extend(['down']*expr.shape[1])
            expr = pd.concat((expr, expr), axis=1)
        expr.columns = add_level(expr.columns, updown_idx)
        return expr.reindex(orig_order, axis=1, level=1)

def _diff_to_int_updown_notation(_diff):
        int_diff = pd.DataFrame(0, _diff.index, _diff.columns.unique(1))
        int_diff[_diff['up'].values] = 1
        if 'down' in _diff.columns.unique(0):
            int_diff[_diff['down'].values] = -1
        return int_diff

def add_updown_mean(agg):
    agg_mean = agg.groupby(axis=1, level=1, sort=False).mean()
    agg_mean.columns = add_level(agg_mean.columns, 'mean')
    return pd.concat([agg, agg_mean], axis=1)

def _add_log2_z(expr):
    expr = np.log2(expr +1)
    expr.columns = add_level(expr.columns, 'log2', at=1)

    m = expr.values.mean()
    s = expr.values.std()
    z_expr = expr.apply(lambda c: (c-m) /s)
    z_expr.columns = add_level(z_expr.columns, 'z', 1, True)

    return pd.concat((expr, z_expr), axis=1).reindex(expr.columns.unique(0), 
                                                     axis=1, level=0)

def add_level(index, label, at=0, replace=False, name=''):
    index = pd.DataFrame(index=index)
    if replace:
        index.reset_index(level=at, drop=True, inplace=True)
    index[name] = label
    order = list(range(index.index.nlevels))
    order.insert(at, -1)
    return index.set_index(name, append=True).reorder_levels(order).index

def annotate(index):
    ref = pd.read_pickle(config.DPRE_PATH + 'mm10/mm10_ensembl_v92_ensg.gzip')
    try:
        return pd.Index(ref.reindex(index).name.values)
    except KeyError as e:
        logger.error(e)

def get_ensgs(names):
    ref = pd.read_pickle(config.DPRE_PATH + 'mm10/mm10_ensembl_v92_ensg.gzip')
    try:
        return ref.loc[names].ensg.tolist()
    except KeyError as e:
        logger.error(e)

def align_indices(data, order, axis=1):
    for i in range(len(data)):
        if data[i] is not None:
            data[i] = data[i].reindex(order, axis=axis)
    return data

# simplify
def filter_trgs(d, max_n, val_th, single_mean=False, ascending=False):
    keep = d.xs('mean', level=-1).iloc[0]
    keep = keep.sort_values(ascending=ascending)[:max_n]
    if val_th:
        keep = keep[keep /keep[0] >val_th]
    if not single_mean:
        d.loc(1)[~d.columns.isin(keep.index)] = np.nan
    else:
        d = keep
    return d





def plot_max_possible_bar(axes, data, ctrl_lbl, bar_args, 
                        draw_colorbar=False, cb_lbl=None, fig=None, 
                        w=None, h=None):
    axes[0].tick_params(labelleft=True)
    axes[0].set_ylim(0, 1)
    axes[0].set_yticks((.5,))
    axes[0].set_yticklabels((ctrl_lbl,), x=.5)
    ax = axes[1]
    [s.set_visible(True) for s in ax.spines.values()]
    im = ax.imshow(data.values, aspect='auto', **bar_args)
    
    if draw_colorbar:
        at = (config.CB_LEFT_SEC/w, 1 - config.CB_TOP/h, 
                config.CB_WIDTH/w, config.CB_HEIGHT/h)

        cb = ax.figure.colorbar(im, cax=fig.add_axes(at), alpha =.3,
                                orientation='horizontal')
        bar_ticks = (bar_args['vmin'], bar_args['vmax'])
        cb.set_ticks(bar_ticks)
        cb.ax.set_xticklabels(bar_ticks)
        cb.ax.set_xlabel(cb_lbl)
        cb.ax.get_xaxis().set_label_position('top')

def _heatmap_cluster(dat, where, ax, which):
    ax.set_visible(True)
    d = dat.T if which == 'columns' else dat
    Y = pdist(d, metric='euclidean')
    Z = linkage(Y, method='complete', metric='euclidean')
    return dendrogram(Z,
                        count_sort = True,
                        no_labels = True,
                        orientation = where, 
                        labels = d.index, 
                        above_threshold_color = config.dendrogram_colors[0],
                        ax = ax)['ivl']

def _init_figure(fig_widths, fig_heights, nplts, spacers):
    width, height = sum(fig_widths), sum(fig_heights)
    ratio = {'width_ratios': list(map(lambda w: w/width, 
                                        fig_widths[1:-2])),
                'height_ratios': list(map(lambda h: h/height, 
                                        fig_heights[1:-2]))}
    # init figure
    fig, axes = plt.subplots(*nplts, figsize=(width, height), 
                                gridspec_kw=ratio)
    axes = clean_axes(axes)
    
    wspace_prop = spacers[0] /np.array(fig_widths[1:-2]).mean()
    hspace_prop = spacers[1] /np.array(fig_heights[1:-2]).mean()
    adj_args = {'left': fig_widths[0] /width, 
                'wspace': wspace_prop,
                'right': 1 - fig_widths[-1] /width,
                'top': 1 - fig_heights[0] /height, 
                'hspace': hspace_prop,
                'bottom': fig_heights[-1] /height}
    fig.subplots_adjust(**adj_args)

    return fig, axes


def open_file(filename):
    if not (filename.endswith('.png') or filename.endswith('.pdf')):
        filename += config.SAVE_FORMAT
    if filename.endswith('.pdf'):
        return PdfPages(filename)

def save_file(fig, filename=None, pp=None, last=False):
    if pp:
        fig.savefig(pp, format='pdf')
    elif filename:
        fig.savefig(filename, format='png')
    logger.info('Plot saved at {}/{}'.format(os.path.abspath(os.curdir), 
                                             filename))

def clean_axes(axes):
    for ax in axes.flatten():
        [s.set_visible(False) for s in ax.spines.values()]
        ax.tick_params(bottom=False, left=False, labelbottom=False, 
                       labelleft=False)
    return axes











# def drop_levels(elements, drop_levels):
#     print(type(elements))
#     for e in elements:
#         for drop in drop_levels:
#             print(e)
#             if type(e) is pd.Series:

#                 e.index = e.index.droplevel(drop)
#             elif type(e) is pd.DataFrame:
#                 e.columns = e.columns.droplevel(drop)
#     return elements


# def norm_proportional_values(result, propof, log=True):
#     propof.index = result.index
#     PLVB = config.PROP_LOW_VALUE_BAHVIOUR
#     toolow = (abs(propof) <PLVB['propof_cap'])
#     if PLVB['b'] == 'cap_result':
#         toolow = (result <PLVB['res_cap']) & toolow.values
        
    # log
    # if log and any(toolow): 
    #     tl_vals = result.mask(~toolow.values)
    #     tl_propof = propof.mask(~toolow.values)
    #     tl_vals = tl_vals.dropna(how='all').dropna(how='all', axis=1)
    #     tl_propof = tl_propof.dropna(how='all').dropna(how='all', axis=1)
    #     f = lambda r, d: ', '.join(r.dropna().round(d).astype(str).tolist())
    #     tl_ns = lambda r: '({}) {}'.format(r.notna().sum(),
    #                                ', '.join(r.dropna().index.tolist()))

    #     tl_res = tl_vals.apply(f, d=1, axis=1)
    #     tl_propof = tl_propof.apply(f, d=3, axis=1)
    #     tl_trgs = tl_vals.apply(tl_ns, axis=1)
    #     df = pd.concat((tl_trgs, tl_propof, tl_res), axis=1)
    #     df.columns = ('targets', 'proportion of', 'result proportion')
    #     tl_n = toolow.sum().sum()
    #     if PLVB['b'] == 'drop':
    #         b = 'drop `result proportion` values ({})'.format(tl_n)
    #     elif PLVB['b'] == 'cap_result':
    #         b = ('substitute `result proportion` values ({}) with defined '
    #             'cap ({})'.format(tl_n, PLVB['res_cap']))
    #     elif PLVB['b'] == 'keep':
    #         b = ('`result proportion` values ({}) will not be droped or '
    #             'substituted'.format(tl_n))
    #     msg = ('Very low proportional values < defined threshold ({}) '
    #            'encountered:\n{}\nSet bahaviour: {}. You can change this '
    #            'behaviour in `config.PROP_LOW_VALUE_BAHVIOUR`.\n'
    #            .format(PLVB['propof_cap'], df.to_string(), b))
    #     logging.warning(msg)
    
    # if PLVB['b'] == 'cap_result':
    #     result[toolow] = PLVB['res_cap']
    # elif PLVB['b'] == 'drop':
    #     result[toolow] = np.nan
    # return result

    
# def _standardized(comb):
#     comb_sdd = comb.copy()
#     comb_sdd[comb_sdd > 1] = 1
#     comb_sdd[comb_sdd < -1] = -1
#     return comb_sdd


# def _update_args(kwargs, newargs):
#     c_kwargs = kwargs.copy()
#     for n_arg, n_arg_val in newargs.items():
#         if n_arg not in kwargs:
#             c_kwargs.update({n_arg: n_arg_val})
#     return c_kwargs

# def norm_to_one(ratio):
#     ratio_normed = [r /sum(ratio) for r in ratio]
#     return ratio_normed
