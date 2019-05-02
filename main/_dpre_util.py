import pandas as pd
import numpy as np
import os
import sys
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from matplotlib.patches import Patch

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

def _diff_to_int_updown_notation(_diff, merge_updown=True):
        int_diff = _diff.astype(int)
        if 'down' in _diff.columns.unique(0):
            int_diff['down'] *= -1
        if merge_updown and 'down' in _diff.columns.unique(0):
            int_diff['up'] = int_diff['up'].mask(_diff['down'], -1)
            int_diff['down'] = int_diff['down'].mask(_diff['up'], 1)
        return int_diff

def add_mgtmean(agg):
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
    path = os.path.dirname(__file__)
    ref = pd.read_pickle(path + '/../mm10/mm10_ensembl_v92_ensg.gzip')
    try:
        return pd.Index(ref.reindex(index).name.values)
    except KeyError as e:
        logger.error(e)
        sys.exit(1)


def get_ensgs(names):
    path = os.path.dirname(__file__)
    ref = pd.read_pickle(path + '/../mm10/mm10_ensembl_v92_ensg.gzip')
    return ref.index[ref.name.isin(names)]

def align_indices(data, order, axis=1):
    for i in range(len(data)):
        if data[i] is not None:
            data[i] = data[i].reindex(order, axis=axis)
    return data

def _make_title(differential, proportional, which, el1, el2, pref='', postf=''):
    which_title = '(euclidean)' if which == 'euclid' else '(intersection)'
    if differential:
        dtype = 'change in '
        if proportional:
            dtype = 'proportional ' + dtype
    elif which == 'euclid' and not differential:
        dtype = 'absolute '
        which_title = ''
    title = '{}{}{}transcriptional similarity {}\nof {} & {}'
    return title.format(pref, dtype, postf, which_title, el1, el2)



    



def plot_required_effect_bar(axes, data, ctrl_lbl, bar_args, 
                        draw_colorbar=False, cb_lbl=None, fig=None, pivot=None, 
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
        if pivot:
            cb.ax.tick_params(labelrotation=90)
        cb.ax.set_xlabel(cb_lbl)
        cb.ax.get_xaxis().set_label_position('top')

def _heatmap_cluster(dat, where, ax, which):
    ax.set_visible(True)
    d = dat.T if which == 'columns' else dat
    Y = pdist(d, metric='euclidean')
    Z = linkage(Y, method='complete', metric='euclidean')
    order = dendrogram(Z,
                        count_sort = True,
                        no_labels = True,
                        orientation = where, 
                        labels = d.index, 
                        above_threshold_color = config.dendrogram_colors[0],
                        ax = ax)['ivl']
    if which == 'rows':
        # for some reason reversed?
        order = order[::-1]
    return order

def _init_figure(fig_widths, fig_heights, nplts, spacers):
    width, height = sum(fig_widths), sum(fig_heights)
    ratio = {'width_ratios': list(map(lambda w: w/width, 
                                        fig_widths[1:-2])),
             'height_ratios': list(map(lambda h: h/height, 
                                        fig_heights[1:-2]))}
    # init figure
    fig, axes = plt.subplots(*nplts, figsize=(width, height), 
                                gridspec_kw=ratio)
    if not isinstance(axes, np.ndarray):
        axes = np.array([axes])
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

def setup_heatmap_xy(x_y, ax, lbls, pivot, hide_lbls, trg_lbl_size, colors):
    
    dim = len(lbls)
    if x_y == 'x':
        # X-axis setup, colorbar bottom
        ax.set_xlim(0, dim)
        ticks = np.arange(.5, dim)
        ax.set_xticks(ticks)
        if not hide_lbls:
            ax.tick_params(labelbottom=True)
            fs = trg_lbl_size*config.FONTS if trg_lbl_size else config.FONTS
            if not pivot:
                ax.set_xticklabels(lbls, rotation=45, ha='right', fontsize=fs, 
                                   rotation_mode='anchor', y=.5)
            else:
                ax.set_xticklabels(lbls, rotation=90, ha='right', va='center', 
                                   fontsize=fs, rotation_mode='anchor', y=.5)

        if colors:
            ax.bar(ticks, 1, 1, color=colors)
    
    elif x_y == 'y':
        ax.set_ylim((-.1, dim +.01))
        ax.set_yticks(np.arange(.5, dim))
        if not hide_lbls:
            ax.tick_params(labelleft=True) 
            if not pivot:
                ax.set_yticklabels(lbls, x=.5)
            else:
                ax.set_yticklabels(lbls, rotation=45, ha='right', x=.5,
                                   rotation_mode='anchor')
        if colors:
            ax.bar(0, 1, color=colors, bottom=np.arange(len(lbls)))

def open_file(filename):
    if not (filename.endswith('.png') or filename.endswith('.pdf')):
        filename += '.'+config.SAVE_FORMAT
    if filename.endswith('.pdf'):
        return filename, PdfPages(filename)
    else:
        return filename, None

def save_file(fig, filename=None, pp=None):
    if pp:
        fig.savefig(pp, format='pdf')
    elif filename:
        fig.savefig(filename, format='png')

def clean_axes(axes):
    np.array([axes])
    for ax in axes.flatten():
        [s.set_visible(False) for s in ax.spines.values()]
        ax.tick_params(bottom=False, left=False, labelbottom=False, 
                       labelleft=False)
    return axes

def color_legend(colors, labels, filename='color_legend.png'):
    assert len(colors) == len(labels), 'colors and labels differ in length'
    spacer.info('\n\n')
    logger.info('Creasting a color legend...')
    fig, ax = plt.subplots(1, 1, figsize=(2, 2))
    clean_axes(np.array([ax]))
    ax.legend(handles=[Patch(color=colors[i], label=labels[i]) 
                for i in range(len(colors))], loc='center')
    fig.savefig(filename)
    logger.info('Legend saved at {}/{}\n\n'
                .format(os.path.abspath(os.curdir), filename))








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

def check_args(trg, smp, which, differential, proportional, 
               hide_distance_bar=None, reorder_to_distance_bar=None,
               cluster_hmx=None, display_similarity=None):
    def check_which(which, trg, smp, diff):
        if which is None:
            if trg._has_expr and smp._has_expr:
                which = 'euclid'
            elif trg._has_diff and smp._has_diff:
                which = 'intersect'
            else:
                logger.error('Either initiate targets and samples with '
                             'expression or with markergenes and diff genes.')
                sys.exit(1)

        msg = 'The {} were initiated without {} data. Cannot use `{}` similarity.'
        if which not in ('euclid', 'intersect'):
            logger.error('Invalid `which` input: `{}`. Valid are `euclid` and '
                        '`intersect`'.format(which))
        elif (which == 'euclid') and not trg._has_expr:
            logger.error(msg.format('targets', 'expression', 'euclid'))
        elif (which == 'euclid') and not smp._has_expr:
            logger.error(msg.format('samples', 'expression', 'euclid'))
        elif (which == 'intersect') and not trg._has_diff:
            logger.error(msg.format('targets', 'merker gene', 'intersect'))
        elif (which == 'intersect') and not smp._has_diff:
            logger.error(msg.format('samples', 'diff genes', 'intersect'))
        elif which == 'euclid' and diff and not smp._ctrl:
            logger.error('To plot the changes in transcriptional similarity '
                         'with which = `euclid`, the samples must be initiated '
                         'with a control.')
        else:
            return which
        sys.exit(1)

    # checks for all plots
    which = check_which(which, trg, smp, differential)
    if which == 'intersect' and not differential:
        differential = True
        logger.warning('For the `intersect` similarity metric, '
                       'differential cannot be False. Was set to True.')
    if proportional and not differential:
                    proportional = False
                    logger.warning('`proportional` can only be used if '
                                   '`differential` is True aswell. Set to False.')

    # checks for 2 heatmaps
    if which == 'euclid' and not hide_distance_bar and not smp._ctrl:
        hide_distance_bar = True
        logger.warning('`hide_distance_bar` cannot be False '
                    'for which = `euclid` if the samples data is '
                    'initialized without a control. Set to True.')
    if reorder_to_distance_bar and hide_distance_bar:
        reorder_to_distance_bar = False
        logger.warning('When `reorder_to_distance_bar` is True, '
                        '`hide_distance_bar` cannot be True. Set '
                        'to False.')
    if reorder_to_distance_bar and cluster_hmx:
        cluster_hmx = False
        logger.warning('Both `reorder_to_distance_bar` and '
                        '`cluster_genes` were set as True. '
                        '`cluster_genes` will be ignored.')                    

    # checks for target_sim and ranked_sim plots
    val = ['mgs mean', 'mgs up', 'mgs down']
    if display_similarity not in val:
        logger.warning('Invalid input for display_similarity: `{}`. ' 
                           'Valid are {}. Set to default `{}`'
                           .format(display_similarity, val, val[0]))
        display_similarity = val[0] 

    if display_similarity == val[2] and not trg._down_mgs: 
        logger.error('Cannot display down markergene similarity because'
                     ' the targets were not initiated with down '
                     'markergenes.')
        sys.exit(1)
    display_similarity = display_similarity[4:]

    return which, differential, proportional, hide_distance_bar, \
           reorder_to_distance_bar, cluster_hmx, display_similarity

