""" Utility module with various helper functions and subplot generaters"""
import pandas as pd
import numpy as np
import os
import sys
import re 
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import is_color_like
from scipy.cluster.hierarchy import distance
from scipy.cluster.hierarchy import linkage 
from scipy.cluster.hierarchy import dendrogram 
from scipy.cluster.hierarchy import optimal_leaf_ordering
from scipy.spatial.distance import pdist

import DPre.main.config as config
from DPre.main._logger import logger, spacer


def _add_mg_types(data, down):
    """Add markergene type index (up and down) to columns at level 0"""
    orig_order = data.columns.unique(0)
    updown_idx = ['up']*data.shape[1]
    if down:
        updown_idx.extend(['down']*data.shape[1])
        data = pd.concat((data, data), axis=1)
    data.columns = _add_level(data.columns, updown_idx)
    return data.reindex(orig_order, axis=1, level=1)

def _bool_to_int_genes(_diff, trans_updown=True, return_merged=False):
    """Take _diff input and convert up-genes to +1, down-genes to -1, optionally 
       transfer up- and down values to each other makeing up- and down subframes
       equal"""
    int_diff = _diff.astype(int)
    if 'down' in _diff.columns.unique(0):
        int_diff['down'] *= -1
    if trans_updown and 'down' in _diff.columns.unique(0):
        int_diff['up'] = int_diff['up'].mask(_diff['down'], -1)
        int_diff['down'] = int_diff['down'].mask(_diff['up'], 1)
    if not return_merged:
        return int_diff
    else:
        return int_diff.xs('up', 1, 0)

def _add_mgtmean(agg):
    """Prodcue the mean between aggregated up- and down mg similarity values"""
    agg_mean = agg.groupby(axis=1, level=1, sort=False).mean()
    agg_mean.columns = _add_level(agg_mean.columns, 'mean')
    return pd.concat([agg, agg_mean], axis=1)

def _add_log2_z(expr, rowwise_sd=False):
    """Compute log2 and z-transformed expression data. Substitute read count in
       expression. Optionally, compute the standad deviation row(gene)-wise. 
       Is used for large datasets like reference transcriptome libraries. 
    """
    expr = np.log2(expr +1)
    expr.columns = _add_level(expr.columns, 'log2', at=1)

    m = expr.values.mean()
    s = expr.values.std() if not rowwise_sd else expr.std(1)
        
    z_expr = expr.apply(lambda c: (c-m) /s)
    z_expr.columns = _add_level(z_expr.columns, 'z', 1)
    return pd.concat((expr, z_expr), axis=1).reindex(expr.columns.unique(0), 
                                                     axis=1, level=0)

def _add_level(index, label, at=0, replace=False, name=''):
    """Add a level with labels 'label' to a pd.MultiIndex"""
    index = pd.DataFrame(index=index)
    if replace:
        index.reset_index(level=at, drop=True, inplace=True)
    index[name] = label
    order = list(range(index.index.nlevels))
    order.insert(at, -1)
    return index.set_index(name, append=True).reorder_levels(order).index

def _get_gene_ann(species):
    """Open the gene annotation reference file (mouse/ human) and return it"""
    path = os.path.dirname(__file__)
    if species == 'mouse':
        return pd.read_pickle(path + '/../gene_ann/mg_ensembl96_GRCm38.p6.gzip')
    elif species == 'human':
        return pd.read_pickle(path + '/../gene_ann/hg_GRCh38.p12.gzip')
    else:
        logger.info('')
        logger.error('Invalid input for species: `{}`. Valid are `mouse` and '
                     '`human`'.format(species))
        sys.exit(1)

def annotate(ensgs, species):
    """ Annotate mouse or human ensg keys. Return the gene names.

    DPre references the ensembl gene annotation v.96 located at 
    DPre/gene_ann. 

    Args:
        ensgs (list, pandas.Index): The collection of ensg keys to annotate
        species (str): The origin species of the genes, 'mouse' or 'human'.
    
    Returns:
        annotated pandas.Index
    """
    ref = _get_gene_ann(species)
    try:
        return pd.Index(ref.reindex(ensgs).name.values)
    except Exception as e:
        logger.error('{}\nDPre references the ensembl gene annotaiton v.96. '
                     'Differently annotated datasets may cause problems.'
                     .format(e))
        sys.exit(1)

def get_ensgs(names, species):
    """ Return the ensg keys for a list of gene names.

    DPre references the ensembl gene annotation v.96 located at 
    DPre/gene_ann. If a gene name has multiple ensg keys, this gene will appear 
    last in the DataFrame regardless of the input order.
    
    Args:
        names (list pandas.Index): The collection of names to return ensg keys 
            for
        species (str): The origin species of the genes, 'mouse' or 'human'.
    
    Returns:
        pandas.Index of ensg keys 
    """
    ref = _get_gene_ann(species)
    try:
        ann = ref.reindex(ref.index[ref.name.isin(names)]).reset_index()
        if ann.name.duplicated().any():
            dupl = pd.Index(ann.name).duplicated()
            ann_dr = ann[~dupl]
            ann_du = ann[dupl]
            ann_dr = ann_dr.set_index('name').reindex(names).reset_index()
            ann_dr.rename({'index': 'name'}, axis=1, inplace=1)
            ann = ann_dr.append(ann_du, sort=False)
            ann.index = np.arange(ann.shape[0])
        else:
            ann = ann.set_index('name').reindex(names).reset_index()
            ann.rename({'index': 'name'}, axis=1, inplace=1)
        return ann
    except Exception as e:
        logger.error('{}\nDPre references the ensembl gene annotaiton v.96. '
                     'Differently annotated datasets may cause problems.'
                     .format(e))
        sys.exit(1)

def _align_indices(data, order, axis=1):
    """Align the indices/ columns in a collection of pandas objects to order"""
    for i in range(len(data)):
        if data[i] is not None:
            data[i] = data[i].reindex(order, axis=axis)
    return data


def _init_figure(fig_widths, fig_heights, nplts, spacers):
    """Calculate the size proportion of each plot element, create figure"""
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
    axes = _clean_axes(axes)
    
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

def _open_file(filename):
    """Open a file based on the filename ending or if not present
       on config.SAVE_FORMAT"""
    if not '.' in filename:
        filename += '.' + config.SAVE_FORMAT
    if filename.endswith('.pdf'):
        return filename, PdfPages(filename)
    else:
        return filename, None
def _save_file(fig, filename=None, pp=None):
    """Save pdf if pp is passed, otherwise use filename to save a .png"""
    if pp:
        fig.savefig(pp, format='pdf')
    elif filename:
        replace = ['$\\mathit{', '}$']
        for repl in replace:
            filename = filename.replace(repl, '')
        fig.savefig(filename)
        plt.close(fig)

def _clean_axes(axes):
    """Remove all spines, ticks and tickalabels"""
    np.array([axes])
    for ax in axes.flatten():
        [s.set_visible(False) for s in ax.spines.values()]
        ax.tick_params(bottom=False, left=False, labelbottom=False, 
                       labelleft=False)
    return axes

def _make_title(differential, metric, el1, el2, pref='', postf=''):
    """Produce the plot title based on plot parmaters, pref and posf are used
       for plot specific adjustments; return the title string"""

    metric_title = 'metric: '
    if metric == 'euclid':
        metric_title += 'L1 Euclidean distance'
    elif metric == 'pearson':
        metric_title += 'Perason correlation' 
    elif metric == 'cosine':
        metric_title += 'cosine distance' 
    elif metric == 'intersect':
         metric_title += 'marker gene intersect'

    if differential:
        dtype = 'Change in '
    else:
        dtype = 'Absolute '
    if pref:
        dtype = dtype.lower()
    title = ('{}{}{}transcriptional similarity \nof {} & {}\n{}'
              .format(pref, dtype, postf, el1, el2, metric_title))
    return title[0].upper() + title[1:]

def _heatmap_cluster(dat, where, ax, metric):
    """Cluster the columns or index with scipy; return the new order"""
    ax.set_visible(True)
    d = dat.T if metric == 'columns' else dat 
    Y = pdist(d, metric='euclidean')
    Z = linkage(Y, method='complete', metric='euclidean')
    order = dendrogram(Z,
                        count_sort = True,
                        no_labels = True,
                        orientation = where, 
                        labels = d.index, 
                        above_threshold_color = config.dendrogram_colors[0],
                        ax = ax)['ivl']
    if metric == 'rows':
        # for some reason reversed?
        order = order[::-1]
    return order

def _plot_distance_bar(axes, data, ctrl_lbl, bar_args, draw_colorbar=False, 
                       cb_lbl=None, fig=None, pivot=None, w=None, h=None):
    """Draw the distance bar on top of the heatmap"""
    # set ylabel on the left 
    axes[0].tick_params(labelleft=True)
    axes[0].set_ylim(0, 1)
    axes[0].set_yticks((.5,))
    axes[0].set_yticklabels((ctrl_lbl,), x=.5)
    # draw the heatmap
    ax = axes[1]
    [s.set_visible(True) for s in ax.spines.values()]
    im = ax.imshow(data.values, aspect='auto', **bar_args)

    # setup the colorbar legend    
    if draw_colorbar:
        at = (config.CB_LEFT_SEC/w, 1- config.CB_TOP/h, config.CB_WIDTH/w, 
              config.CB_HEIGHT/h)
        cb = ax.figure.colorbar(im, cax=fig.add_axes(at), alpha =.3,
                                orientation='horizontal')
        bar_ticks = (bar_args['vmin'], bar_args['vmax'])
        cb.set_ticks(bar_ticks)
        cb.ax.set_xticklabels(bar_ticks)
        if pivot:
            cb.ax.tick_params(labelrotation=90)
        cb.ax.set_xlabel(cb_lbl)
        cb.ax.get_xaxis().set_label_position('top')

def _setup_heatmap_xy(x_y, ax, lbls, pivot, hide_lbls, lbl_size, colors):
    """Setting all paramters for the x- and y axis of the two heatmap plots"""
    dim = len(lbls)
    if x_y == 'x':
        # X-axis setup, colorbar bottom
        ax.set_xlim(0, dim)
        ticks = np.arange(.5, dim)
        ax.set_xticks(ticks)
        if not hide_lbls:
            ax.tick_params(labelbottom=True)
            fs = lbl_size*config.FONTS if lbl_size else config.FONTS
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
            fs = lbl_size*config.FONTS if lbl_size else config.FONTS
            if not pivot:
                ax.set_yticklabels(lbls, x=.5, fontsize=fs)
            else:
                ax.set_yticklabels(lbls, rotation=45, ha='right', x=.5,
                                   fontsize=fs, rotation_mode='anchor')
        if colors:
            ax.bar(0, 1, color=colors, bottom=np.arange(len(lbls)))

def _check_args(trg, smp, metric, differential, 
                hide_distance_bar=None, reorder_to_distance_bar=None,
                cluster_hmx=None, display_similarity=False):
    """General purpose plot argument checker; returns (modified) input values"""
    def check_metric(metric, trg, smp, diff):
        # check if the samples and targets have equivalent data to compare
        if metric is None:
            if trg._has_expr and smp._has_expr:
                metric = 'euclid'
            elif trg._has_diff and smp._has_diff:
                metric = 'intersect'
            else:
                logger.error('Either initiate targets and samples with '
                             'expression or with markergenes and diff genes.')
                sys.exit(1)
        msg = 'The {} were initiated without {} data. Cannot use `{}` similarity.'
        if metric not in ('euclid', 'intersect', 'cosine', 'pearson'):
            logger.error('Invalid `metric` input: `{}`. Valid are `euclid` and '
                        '`intersect`'.format(metric))
        elif (metric == 'euclid') and not trg._has_expr:
            logger.error(msg.format('targets', 'expression', 'euclid'))
        elif (metric == 'euclid') and not smp._has_expr:
            logger.error(msg.format('samples', 'expression', 'euclid'))
        elif (metric == 'intersect') and not trg._has_diff:
            logger.error(msg.format('targets', 'merker gene', 'intersect'))
        elif (metric == 'intersect') and not smp._has_diff:
            logger.error(msg.format('samples', 'diff genes', 'intersect'))
        elif metric == 'euclid' and diff and not smp._ctrl:
            logger.error('To plot the changes in transcriptional similarity '
                         'with metric = `euclid`, the samples must be initiated '
                         'with a control.')
        else:
            # valid metric
            return metric
        # invalid metric
        sys.exit(1)

    # checks for all plots
    metric = check_metric(metric, trg, smp, differential)
    if metric == 'intersect' and not differential:
        differential = True
        logger.warning('For the `intersect` similarity metric, '
                       'differential cannot be False. Was set to True.')
    # if proportional and not differential:
    #                 proportional = False
    #                 logger.warning('`proportional` can only be used if '
    #                                '`differential` is True aswell. Set to False.')

    # checks for 2 heatmaps
    if metric == 'euclid' and not hide_distance_bar and not smp._ctrl:
        hide_distance_bar = True
        logger.warning('`hide_distance_bar` cannot be False '
                    'for metric = `euclid` if the samples data is '
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

    if display_similarity is not False:
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

    return metric, differential, hide_distance_bar, reorder_to_distance_bar, \
           cluster_hmx, display_similarity

def plot_color_legend(labels, colors, ncolumns=4, filename='color_legend.png'):
    """Plot a custom color legend.
    
       Takes a list of labels and colors and links them to produce a color 
       legend. Useful for marking  sub-groups in samples/ targets elements.

    Args:
        labels (list): the list of labels in the legend
        colors (list): the list of colors correspoding to the labels. Colors 
            must be interpretable by matplotlib: for example, 'w', #ffffff, 
            (1,1,1) all refer to white.
        filename (str, optional): the filename to save the legend. Defaults to
            './color_legend.png'
        ncolumns (int, optional): the number of columns in the legend. Defaults 
            to 1.
    """
    spacer.info('\n\n')
    assert len(colors) == len(labels), 'colors and labels differ in length'
    inv_cols = [c for c in colors if not is_color_like(c)]
    if inv_cols:
        logger.error('The following colors are not recognized as colors by '
                     'matplotlib: {}'.format(inv_cols))
        sys.exit(1)
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    _clean_axes(np.array([ax]))
    ax.legend(handles=[Patch(color=colors[i], label=labels[i]) 
                for i in range(len(colors))], loc='center', ncol=ncolumns)
    fig.savefig(filename)
    plt.close()
    logger.info('Color legend generated and saved at {}/{}'
                .format(os.path.abspath(os.curdir), filename))