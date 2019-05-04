import pandas as pd
import numpy as np
import sys
import os
import copy
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D

from DPre.main._differential import _differential
from DPre.main.samples import Samples
from DPre.main._logger import spacer, logger, log_plot
import DPre.main.config as config
import DPre.main._dpre_util as util

class Targets(_differential):
    """class describing the data to compare similarity against

    Targets can hold lists of markergenes and expression data identifying a 
    collection of comparitive transcriptional identities, the targets. 
    Targets are recomended to be initiated with markergenes. 

    Arguments:
        markergenes (optional): directory with deseq2 output OR directories with 
            up- and (optional) down genelists OR pandas.DataFrame. Defaults to None.
            Genelist data have an ensg key in the first column or contain an 'ensg' 
            column label.
            When passing a DataFrame, the index should consist of all genes, 
            the columns of a pandas.MultiIndex with up or up and down at level 0 and 
            the element names at level 1. THe dtype is bool, markergenes are marked 
            as True values. 
        expression (optional): directory with expression tsv file OR 
            pandas.DataFrame. Defaults to None.
            The .tsv input should have an ensg key in the first column or ensg 
            column label. Columns `loc`, `name`, `tss_loc` and `strand` are removed
            automatically. The data should be exclusively numerical without NaN's.
            When passing a DataFrame, the data should be log2- and z-transformed. 
            The index should consist of all genes, the columns of a 
            pandas.MultiIndex with the element names at level 0, and `log2` & `z` at 
            level 1.
        ignore_down_mgs (bool, optional): Even if found in markergene input, do not
            use the down markergenes for the analysis. Defaults to False.
        override_namematcher (bool, optional): When both markergenes and expression 
            passed, this overrides the element names in markergenes. Defaults to 
            False. When False, element names in markergenes and expression are 
            expected to match perfectly.
        name (str, optional): Name label of the Targets. Defaults to 'Targets'. Used
            in logging and plot annotations.
        log: (bool, optional): Log the Targets initiation. Defaults to True.
        
    Note:
        At least one of markergenes and expression must be passed. When both are 
        passed, the inputs must have the same elment order. Markergene input is
        alphabetically sorted. It is recommended to initiate with markergenes.
    """
    def __init__(self, markergenes=None, expression=None, name=None, 
                 ignore_down_mgs=False, override_namematcher=False, 
                 species='mouse', log=True):
        # call _differential __init__ method
        super().__init__(diff_genes=markergenes, expression=expression, 
                         name=name, override_namematcher=override_namematcher, 
                         log=log)
        # define if down markergenes are used, which species and the overlap dict
        self._down_mgs = not ignore_down_mgs
        self._species = species
        self._overlaps = {}
        
        # remove down mgs from _diff if in there but not desired by the user
        if self._has_diff:
            if not self._down_mgs and 'down' in self._diff.columns.unique(0):
                self._diff.drop('down', axis=1, level=0, inplace=True)
            elif 'down' not in self._diff.columns.unique(0) and self._down_mgs:
                self._down_mgs = False
            if log:
                spacer.info('')
                n = self._diff.sum().unstack(0).reindex(self.names).to_string()
                logger.info('Number of markergenes: \n{}'.format(n))
        # inform that not passing markergenes is not recommended
        else:
            self._down_mgs = False
            spacer.warning('')
            logger.warning('The targets `{}` are initiated without markergenes.'
                           ' Note that comparing against all genes can '
                           'lead to low accuracy for defining transcriptional '
                           'similarity.'.format(self.name))
        
        # store a masked version of _expr that holds only the markergenes
        if self._has_expr:
            expr_mgs = util._add_mg_types(self._expr.copy(), self._down_mgs)
            if self._has_diff: 
                mg_mask = lambda trg: trg.mask(~self._diff[trg.columns[0][:-1]])
                expr_mgs = expr_mgs.groupby(level=(0,1), axis=1).apply(mg_mask)
                self._expr_mgs = expr_mgs.reindex(self._mgs)
            else:
                self._expr_mgs = expr_mgs
        spacer.info('\n\n')
        self._log_init(log)

    def __repr__(self):
        return ('\n=|=|= Targets-instance =|=|=\nname = {};\nelements = {};\n'
                'n = {};\nmarkergenes data = {};\nexpression data = {}\n'
                .format(self.name, self.names, len(self), self._has_diff, 
                        self._has_expr))
    @property
    def _mgs(self):
        if self._has_diff:
            return self._diff[self._diff.any(1)].index
        elif self._has_expr:
            return self._expr[self._expr.any(1)].index

    @property
    def _mg_types(self):
        return ['up', 'down'] if self._down_mgs else ['up']

    def _overlap_samples(self, samples, which):
        spacer.info('')
        logger.info('Overlappping samples ({}): `{}` on Targets: `{}` ... '
                    .format(which, samples.name, self.name))
        
        # check marker gene detection before overlapping
        det = self.detec_mgs_prop(samples, filename=None).reindex(self.names, 
                                                                  level=1)
        keep = det[det.proportion > config.DROP_TARGET_DETEC_THR].index.unique(1)
        if len(keep) != len(self):
            # drop targets with too few detected genes
            dr = pd.Index(self.names).difference(keep).tolist()
            logger.info('{} target elements dropped due to marker gene detection'
                        ' proportions lower {} (Set in config.DROP_TARGET_'
                        'DETEC_THR):\n{}'
                        .format(len(dr), config.DROP_TARGET_DETEC_THR, dr))
            self = self.slice_elements(keep, log=False)

        if which == 'euclid':
            # get the z expression of target markergenes and samples (all)
            # if the targets have down markergenes, duplicate sample data
            trg_data = self._expr_mgs.xs('z', 1, 2)
            smp_data = samples._expr.xs('z', 1, 1)
            smp_data = util._add_mg_types(smp_data, self._down_mgs)
            subst = samples._expr.xs('z', level=1, axis=1).min().min()
        elif which == 'intersect':
            # substitute diff data with +1 for upregualted genes, -1 for down
            # if the targets have no down markergenes, slice the sample data
            # same as for samples, but don't merge up and down markergene lists
            smp_data = util._diff_to_int_updown_notation(samples._diff, True)
            smp_data = smp_data[self._mg_types]
            diff_mgs = self._diff.reindex(self._mgs)
            trg_data = util._diff_to_int_updown_notation(diff_mgs, False)
            # non markergenes are set to NaN 
            trg_data.mask(trg_data == 0, inplace=True)
            subst = 0

        # multiply the sample data len(target data) times for groupby operation
        smp_data = smp_data.reindex(self._mgs)
        smp_exts = [d for _, dd in smp_data.iteritems() 
                    for d in [dd]*len(self.names)]
        smp_ext = pd.concat(smp_exts, axis=1)
        lvls = (self._mg_types, smp_ext.columns.unique(1), self.names)
        smp_ext.columns = pd.MultiIndex.from_product(lvls).swaplevel()
        
        # ignore not detected marker genes 
        smp_ext.dropna(inplace=True)
        trg_data = trg_data.reindex(smp_ext.index)

        # core
        # calculate the differnece between the z target marker gene expression 
        # in the targets and samples. positive value = sample higher expressed
        eucl_dist = lambda smp_d: smp_d -trg_data[smp_d.columns.unique(0)].values 
        # calculate matches between targets and samples. A match results in 
        # -1+-1= -2 or 1+1=2, mismatches in 0, no diff-gene/ no marker gene in 1 
        # Matches are set to 1, mismatches to -1, no overlap to 0
        def mg_inters(smp_d):
            m = abs(smp_d + trg_data[smp_d.columns[0][0]].values)
            return m.mask(m == 2, 1).mask(m == 1, 0).mask(m == 0, -1)
        do_ovp = eucl_dist if which == 'euclid' else mg_inters
        ovp = smp_ext.groupby(axis=1, level=(0,2), sort=False).apply(do_ovp)
        self._overlaps['{}-{}'.format(id(samples), which)] = ovp
        return ovp

    def detec_mgs_prop(self, samples, filename='detec_mgs_prop.png'):
        # get proportion of detected markergenes
        if samples._has_expr:
            smp_d = samples._expr.reindex(self._mgs).notna().iloc(1)[0]
        elif samples._has_diff:
            smp_d = samples._diff.reindex(self._mgs).notna().iloc(1)[0]
        if self._has_diff:
            trg_d = self._diff
        elif self._has_expr:
            cols = pd.MultiIndex.from_product((self._mg_types, self.names))
            trg_d = pd.DataFrame(True, self._expr.index, cols)
        det = trg_d.reindex(self._mgs).apply(lambda trg: trg & smp_d).sum()
        n_mgs = trg_d.sum()
        order = (det/n_mgs).sort_values().index
        # log proportion of detected markergenes
        df = pd.DataFrame({'n markergenes': n_mgs.reindex(order), 
                           'detected in samples': det.reindex(order).values, 
                           'proportion': (det/n_mgs).reindex(order).values})
        n_trgs = 10 if not len(order) <20 else int(len(order)/2)
        edges = order.droplevel(0)[:n_trgs].append(order.droplevel(0)[-n_trgs:])
        logger.info('Detection of target ({}) marker genes in sample data ({}):'
                    '\n{}\nShown are the {} edge proportion values.'
                    .format(self.name, samples.name, 
                            df.loc[(slice(None), edges), :].to_string(), 
                            len(edges)))
        if filename:
            fig, ax = plt.subplots()
            ax.bar(np.arange(len(order)), df.proportion, width=1, edgecolor='k',
                              color=self.get_colors(order.get_level_values(1)))
            ax.hlines(config.DROP_TARGET_DETEC_THR, 0, len(self))
            ax.set_xlabel(self.name, fontsize=4)
            ax.set_ylabel('Proportion of detected markergenes', fontsize=4)
            ax.set_title('Proportion of detected {} markergenes in {}\nline = drop'
                      ' threshhold'.format(self.name, samples.name), fontsize=6)
            fig.savefig(fname=filename)
            logger.info('Plot saved at {}\n'
                        .format(os.path.abspath(filename)))
            plt.close()
        return df

    def _get_from_overlap(self, samples, which, genes_diff=False, genes_prop=False,
                          genes_agg_diff=False, genes_agg_prop=False, 
                          drop_ctrl=True, inters_to_updown_not=True, log=True):
        # check if overlap has already been computed, if not do overlap
        try:
            ovp = self._overlaps['{}-{}'.format(id(samples), which)].copy()
        except KeyError:
            ovp = self._overlap_samples(samples, which).copy()

        if log:
            logger.info('Selecting and processing overlap data...')
        t_ord = ovp.columns.unique(1)
        s_ord = ovp.columns.unique(2)
        if which == 'euclid':
            # get single gene control ovp values and aggrevate (absolute mean)
            if samples._ctrl:
                ctrl_genes = ovp.xs(samples._ctrl, 1, 2, False)
                c_a = ctrl_genes.abs().mean().unstack().T
                ctrl_agg = util.add_mgtmean(c_a.reindex(t_ord, axis=1, level=1))
            else:
                ctrl_genes = None
                ctrl_agg = None

            if not genes_diff:
                genes = ovp
            # get samples single gene effects
            else:
                def gene_diff(smp_ovp, alter=None):
                    ctrl = smp_ovp.xs(samples._ctrl, 1, 2).iloc(1)[0]
                    to_diff = lambda smp:  ctrl.abs() - smp.abs()
                    diff = smp_ovp.apply(to_diff)
                    if genes_prop:
                        to_prop = lambda smp: smp / ctrl.abs()
                        # make interpretable, 1 = 100% -1 = -100%
                        prop_interpr = ((diff.apply(to_prop)-1).abs() *-1) +1
                        # cap negative single gene proportional values
                        return prop_interpr.mask(prop_interpr < -3, -3)
                    return diff
                genes = ovp.groupby(level=(0,1), axis=1, sort=False).apply(gene_diff)
            
            # aggrevate the sample data
            if not genes_agg_diff:
                agg = ovp.abs().mean().unstack((0,1)).reindex(s_ord)
                agg = util.add_mgtmean(agg)
            else:
                # get aggreavated absolute mean effects
                def agg_diff(smp_ovp):
                    agg_mean = smp_ovp.abs().mean()
                    ctrl_mean = agg_mean.xs(samples._ctrl, level=2).values
                    if genes_agg_diff:
                        agg_mean = ctrl_mean - agg_mean
                        if genes_agg_prop:
                            agg_mean /= ctrl_mean
                    return agg_mean.unstack((0,1)).reindex(agg_mean.index.unique(2))
                agg = ovp.groupby(level=(0,1), axis=1, sort=False).apply(agg_diff)
                agg = util.add_mgtmean(agg.droplevel((0,1), axis=1))
            
            if samples._ctrl and drop_ctrl:
                genes.drop(samples._ctrl, axis=1, level=2, inplace=True)
                agg.drop(samples._ctrl, inplace=True)
            return genes, agg, ctrl_genes, ctrl_agg
        
        elif which == 'intersect':
            # ovp stores matches (1) and mismatches (-1). This option is for 
            # reverting that notation for down mgs, matches -1, mismatches 1
            n_mgs = ovp.notna().sum().xs(s_ord[0], level=2)
            n_mgs.index = util._add_level(n_mgs.index, 'n_mgs', at=2)
            n_mgs = util.add_mgtmean(n_mgs.unstack((0,1))).astype(int)

            agg = ovp.sum().unstack((0,1)).reindex(s_ord)
            agg = util.add_mgtmean(agg)
            if genes_agg_prop:
                agg /= n_mgs.iloc[0]

            if samples._ctrl and drop_ctrl:
                ovp.drop(samples._ctrl, axis=1, level=2, inplace=True)
                agg.drop(samples._ctrl, inplace=True)
            if inters_to_updown_not and self._down_mgs:
                ovp['down'] *= -1
            return ovp, agg, None, n_mgs

    def target_similarity_heatmap(self, 
                                  # plot data
                                  samples, 
                                  which, 
                                  differential = False,
                                  proportional = False, 
                                  display_similarity = 'mgs mean',
                                  # data ordering
                                  cluster_targets = False,
                                  cluster_samples = False,
                                  reorder_to_distance_bar = False,
                                  # general settings
                                  pivot = False,
                                  heatmap_width = None,
                                  heatmap_height = None,
                                  heatmap_range = None,
                                  distance_bar_range = None,
                                  targetlabels_space = None,
                                  samplelabels_space = None,
                                  targetlabels_size = None,
                                  title = True, 
                                  # show/ hide elements 
                                  hide_colorbar_legend = False,
                                  hide_distance_bar = False,
                                  hide_targetlabels = False,
                                  hide_targets_dendrogram = False,
                                  hide_targets_colorbar = False,
                                  hide_samplelabels = False,
                                  show_samples_dendrogram = False,
                                  show_samples_colorbar = False,
                                  filename = 'target_similarity_hm.png'):
        # check user input for errors and incompatibilities
        def _check_args():
            nonlocal which
            nonlocal differential
            nonlocal proportional

            nonlocal cluster_targets
            nonlocal reorder_to_distance_bar
            nonlocal hide_distance_bar
            nonlocal display_similarity

            # check general basic input requirements
            a = util._check_args(self, samples, which, differential, proportional, 
                                hide_distance_bar, reorder_to_distance_bar,
                                cluster_targets, display_similarity)
            which, differential, proportional, hide_distance_bar, \
            reorder_to_distance_bar, cluster_targets, display_similarity = a
            logger.info('Arguments passed. Getting data now ...')

        # get the specific overlap data, plot the mean of up and down mgs
        def get_data():
            _, sim, _, ctrl_sim = self._get_from_overlap(samples, which, 
                                                    genes_agg_diff=differential, 
                                                    genes_agg_prop=proportional)
            if ctrl_sim is None:
                ctrl_sim =  pd.DataFrame(0, [0], sim.columns)
            return [sim.xs(display_similarity, 1, 0), 
                    ctrl_sim.xs(display_similarity, 1, 0)]

        # get plot lims
        def get_caps():
            mini = abs(data[0].min().min())
            maxi = abs(data[0].max().max())
            cap = round(max((mini, maxi)), 1)
            re_cap = round(data[1].iloc[0].max(), 1)
            if which == 'euclid' and not differential:
                cap = re_cap = max((cap, re_cap))
            return cap, re_cap
        
        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            nplts = [4,3]

            fig_widths = [.0001] *(nplts[1] +3)
            l = config.HM_LEFT if not pivot else config.HM_PIVOT_LEFT
            fig_widths[0] = samplelabels_space if samplelabels_space else l
            if show_samples_colorbar:
                fig_widths[1] = config.HM_Y_COLORBAR
            fig_widths[2] = config.HM_SQUARE_SIZE * len(self)
            if heatmap_width:
                fig_widths[2] *= heatmap_width
            if cluster_samples and show_samples_dendrogram:
                fig_widths[3] = config.HM_Y_DENDROGRAM
            fig_widths[4] = config.HM_WSPACE * (nplts[1]-1)
            fig_widths[5] = config.HM_RIGHT

            fig_heights = [.0001] *(nplts[0] +3)
            fig_heights[0] = config.HM_TOP
            if cluster_targets and not hide_targets_dendrogram:
                fig_heights[1] = config.HM_X_DENDROGRAM
            if not hide_distance_bar:
                fig_heights[2] = config.HM_REQU_EFF_BAR
            fig_heights[3] = config.HM_SQUARE_SIZE * len(samples._names_noctrl)
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if not hide_targets_colorbar:
                fig_heights[4] = config.HM_X_COLORBAR
            fig_heights[5] = config.HM_HSPACE * (nplts[0]-1)
            b = config.HM_BOTTOM if not pivot else config.HM_PIVOT_BOTTOM
            fig_heights[6] = targetlabels_space if targetlabels_space else b

            return nplts, fig_widths, fig_heights

        # draw plot
        def do_plot():
            width, height = sum(fig_widths), sum(fig_heights)
            fig, axes = util._init_figure(fig_widths, fig_heights, nplts, 
                                       (config.HM_WSPACE, config.HM_HSPACE))
            sim, ctrl_sim = data

            # set plot title
            if title in ('None', 'none', 'False', 'false', 'F', 'f'):
                this_t = False
            if title:
                if title and type(title) is not str:
                    this_t = util._make_title(differential, proportional, which, 
                                              samples.name, self.name)
                else:
                    this_t = title
                if not pivot:
                    fig.suptitle(this_t, fontsize=config.FONTS *1.1,
                                 y=1-((config.HM_TOP/height *.1)))
                else:
                    axes[2, 0].set_ylabel(this_t, fontsize=config.FONTS *1.1,
                                          labelpad=30)


            # cluster targets/ samples and draw dendrograms
            if cluster_targets:
                at = axes[0, 1] if not hide_targets_dendrogram else axes[0, 0]
                order = util._heatmap_cluster(sim, 'top', at, 'columns')
                sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
            if cluster_samples:
                at = axes[2, 2] if show_samples_dendrogram else axes[0, 0]
                order = util._heatmap_cluster(sim, 'right', at, 'rows')
                sim = sim.reindex(order)
            axes[0, 0].set_visible(False)
            
            # draw distance effect bar
            if not hide_distance_bar:
                if reorder_to_distance_bar:
                    order = ctrl_sim.iloc[0].sort_values().index
                    sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
                draw_cb = True if not hide_colorbar_legend else False
                if which == 'euclid' and not differential:
                        draw_cb = False
                ctrl_lbl = ''
                if which == 'euclid' and not hide_samplelabels:
                        ctrl_lbl = samples._ctrl
                bar_args = {'cmap': 'afmhot', 'vmin': 0, 'vmax': re_cap}
                cb_lbl = 'number of markergenes\n' if which == 'intersect' else \
                         'base expr. similarity\n[mean abs. eucl. distance]'
                if distance_bar_range is not None:
                    bar_args.update({'vmin': distance_bar_range[0], 
                                     'vmax': distance_bar_range[1]})
                util._plot_distance_bar(axes[1, :2], ctrl_sim,
                                            ctrl_lbl, bar_args, draw_cb, 
                                            cb_lbl, fig, pivot, width, height)

            # setup heatmap x,y axis, including the colorbars
            cols = self.get_colors(sim.columns) if not hide_targets_colorbar \
                   else None
            util._setup_heatmap_xy('x', axes[3, 1], sim.columns, pivot,
                                  hide_targetlabels, targetlabels_size, cols)
                   
            cols = samples.get_colors(sim.index[::-1]) if show_samples_colorbar \
                   else None
            util._setup_heatmap_xy('y', axes[2, 0], sim.index[::-1], pivot, 
                                  hide_samplelabels, targetlabels_size, cols)

            ax = axes[2, 1]
            ax.set_yticks(np.arange(0, sim.shape[0]))
            ax.set_xticks(np.arange(0, sim.shape[1]))
            hm_args = {'cmap': 'RdBu_r', 'vmin': -cap, 'vmax': cap}
            if which == 'euclid' and differential and not proportional:
                cb_lbl = config.AGG_EUCL_DIFF_NOPROP
            if which == 'euclid' and differential and proportional:
                cb_lbl = config.AGG_EUCL_DIFF_PROP
            elif which == 'euclid' and not differential:
                hm_args = {'cmap': 'afmhot', 'vmin': 0, 'vmax': cap}
                cb_lbl = config.AGG_EUCL_NODIFF
            elif which == 'intersect' and not proportional:
                cb_lbl = config.AGG_INTE_DIFF_NOPROP
            elif which == 'intersect' and proportional:
                cb_lbl = config.AGG_INTE_DIFF_PROP
            if heatmap_range is not None:
                hm_args.update({'vmin': heatmap_range[0], 
                                'vmax': heatmap_range[1]})
            im = ax.imshow(sim.values, aspect='auto', **hm_args)
            
            # setup heatmap colorbar legend and draw
            if not hide_colorbar_legend:
                at = (config.CB_LEFT/width, 1-config.CB_TOP/height, 
                    config.CB_WIDTH/width, config.CB_HEIGHT/height)
                cax = fig.add_axes(at)
                cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
                
                bar_ticks = [hm_args['vmin'], hm_args['vmax']]
                if which == 'intersect' and not proportional:
                    bar_ticks = [int(bar_ticks[0]), int(bar_ticks[1])]
                cb.set_ticks(bar_ticks)
                cb.ax.set_xticklabels(bar_ticks)
                if pivot:
                    cb.ax.tick_params(labelrotation=90)
                cb.ax.set_xlabel(cb_lbl)
                cb.ax.get_xaxis().set_label_position('top')

            return fig, axes, (sim, ctrl_sim)

        spacer.info('\n\n' + log_plot)
        logger.info('Plot: {} & {}'.format(self.name, samples.name))
        _check_args()
        data = get_data()
        cap, re_cap = get_caps()

        nplts, fig_widths, fig_heights = get_plot_sizes()
        spacer.info('')
        logger.info('Drawing...')
        filename, pp = util._open_file(filename)
        fig, axes, data = do_plot()
        # plt.show()
        util._save_file(fig, filename=filename, pp=pp)
        if pp:
            pp.close()
        logger.info('Plot saved at {}/{}\n\n'
                    .format(os.path.abspath(os.curdir), filename))
        return fig, axes, data
        

    def gene_similarity_heatmap(self,
                                # plot data
                                samples,  
                                which,
                                differential = False,
                                proportional = False, 
                                display_genes = 'variant',
                                gene_number = 45,
                                specific_genes = None,
                                custom_target_genelist = None,
                                # data ordering
                                cluster_genes = False,
                                cluster_samples = False,
                                reorder_to_distance_bar = False,
                                # general settings
                                pivot = False,
                                heatmap_width = None,
                                heatmap_height = None,
                                heatmap_range = None,
                                distance_bar_range = None,
                                sum_plot_range = None,
                                genelabels_space = None,
                                genelabels_size = None,
                                samplelabels_space_left = None,
                                samplelabels_space_right = None,
                                title = True, 
                                # show/ hide elements 
                                hide_colorbar_legend = False,
                                hide_distance_bar = False,
                                hide_sum_plot = False,
                                show_samplelabels_sum_plot = False,
                                hide_genelabels = False,
                                hide_genes_dendrogram = False,
                                show_genes_colorbar = None,
                                hide_samplelabels = False,
                                show_samples_dendrogram = False,
                                show_samples_colorbar = False,
                                filename = 'gene_similarity_hm.pdf'):
        
        # check user input for errors and incompatibilities
        def _check_args():
            nonlocal which
            nonlocal differential
            nonlocal proportional
            nonlocal display_genes
            nonlocal specific_genes
            nonlocal custom_target_genelist

            nonlocal hide_distance_bar
            nonlocal reorder_to_distance_bar
            nonlocal cluster_genes
            nonlocal show_genes_colorbar

            # check general basic input requirements
            a = util._check_args(self, samples, which, differential, proportional, 
                                hide_distance_bar, reorder_to_distance_bar,
                                cluster_genes)
            which, differential, proportional, hide_distance_bar, \
            reorder_to_distance_bar, cluster_genes, _ = a

            # check main data input
            if custom_target_genelist is not None and display_genes:
                display_genes = None
                logger.info('Both `display_genes` and '
                            '`custom_target_genelist` were passed. '
                            '`display_genes` will be ignored.')
            if display_genes:
                val = ['variant', 'greatest', 'increasing', 'decreasing']
                if not differential:
                     val = ['variant', 'distant', 'similar']
                if display_genes not in val:
                    logger.error('The passed value for display_genes: `{}` is '
                                 'invalid. Valid options when `differential` is'
                                 ' {} are {}.'
                                 .format(display_genes, differential, val))
                    sys.exit(1)
            elif custom_target_genelist is None and specific_genes is None:
                logger.error('None of `display_genes`, `specific_genes` or '
                             '`custom_target_genelist` were passed')
                sys.exit(1)
            elif custom_target_genelist is not None and specific_genes is not None:
                specific_genes = None
                msg = ('Both `specific_genes` and `custom_target_genelist` were'
                       ' passed. `specific_genes` will be ignored.')
                logger.info(msg)
            
            # modify arguments for convneience 
            if which == 'euclid' and not differential:
                if display_genes == 'distant':
                    display_genes = 'increasing'
                elif display_genes == 'similar':
                    display_genes = 'decreasing'
            if show_genes_colorbar == True:
                if specific_genes:
                    show_genes_colorbar = dict.fromkeys(specific_genes, 
                                                   config.colors[1])
                else:
                    show_genes_colorbar = None

            # get a list of generally valid annotated genes
            genes = pd.DataFrame({'name': util.annotate(self._mgs, self._species), 
                                  'ensg': self._mgs })
            if specific_genes is not None or custom_target_genelist is not None:
                # for gene input check if genes are detected in the target data
                if specific_genes is not None:
                    inp_gl = pd.Index(specific_genes).drop_duplicates()
                    val_gl = pd.Index(genes.name.values)
                    isin = 'markergenes'
                elif custom_target_genelist is not None:
                    inp_gl = pd.Index(custom_target_genelist).drop_duplicates()
                    val_gl = self._detec_genes.intersection(samples._detec_genes)
                    isin = 'detected genes'
                    val_gl = pd.Index(util.annotate(val_gl, self._species))
                
                inv = [g for g in inp_gl if g not in val_gl]
                inp_gl = inp_gl.drop(inv) 
                if inv:
                    logger.warning('{} ({}/{}) are not {} in any of the targets'
                                   ' in `{}`. These genes will not be included.'
                                   .format(inv, len(inv), len(inv)+len(inp_gl), 
                                           isin, self.name))
                    if len(inv) == (len(inv)+len(inp_gl)):
                        sys.exit(1)
                # update passed list
                if specific_genes is not None:
                    specific_genes = inp_gl
                elif custom_target_genelist is not None:
                    genes = pd.DataFrame({'name': inp_gl, 
                                          'ensg': util.get_ensgs(inp_gl, 
                                                                 self._species)})
            logger.info('Arguments passed. Getting data now ...')
            return genes                   

        # get the specific overlap data and pick out the genes to display
        def get_data():
            
            # init a new target where all genes are markergenes of all targets
            if custom_target_genelist:
                nonlocal self
                diff = pd.DataFrame(True, genes.ensg, self.names)
                diff = util._add_mg_types(diff, False)
                args = {'markergenes': diff}
                if which == 'euclid':
                    expr = self._expr.loc[genes.ensg].copy()
                    args.update({'expression': expr})
                self = Targets(name='custom genelist', ignore_down_mgs=True, 
                               log=False, **args)
            sim, _, ctrl_sim, _ = self._get_from_overlap(samples, which, 
                                                      genes_diff=differential, 
                                                      genes_prop=proportional)
            if which == 'euclid':
                if samples._ctrl:
                    ctrl_name = samples._ctrl
                    ctrl_sim = ctrl_sim.abs()
                else:
                    # dummy, won't be used
                    ctrl_name = sim.columns[0][2]
                    ctrl_sim = sim.xs(ctrl_name, 1, 2, False)
                if not differential:
                    sim = sim.abs()
                if proportional:
                    s_noprop, _ , cs_noprop, _ = self._get_from_overlap(samples,
                                                               which, log=False) 

            # init mutable dict with target and markegene type keys
            data = dict((trg, dict((mgt, None) for mgt in self._mg_types))
                        for trg in self.names)
            # select genes, form the 3 data elements target simlarity (heatmap), 
            # ctrl_sim (required_effect_bar), agg_sim (sumplot)
            def sel_genes(trg_sim, genes):
                mgt = trg_sim.columns[0][0]
                trg = trg_sim.columns[0][1]
                get_genes = pd.Index([])
                trg_sim.dropna(inplace=True)
                
                if display_genes:
                    # sort overlap based on passed metric, slice to gene number
                    if display_genes == 'variant':
                        idx = trg_sim.var(1).sort_values().index[::-1]
                    elif display_genes == 'greatest' and which == 'euclid':
                        idx = trg_sim.abs().max(1).sort_values().index[::-1]
                    elif display_genes == 'greatest' and which == 'intersect':
                        idx = trg_sim.abs().sum(1).sort_values().index[::-1]
                    elif display_genes == 'increasing' and (which == 'euclid'):
                        idx = trg_sim.max(1).sort_values().index[::-1]
                    elif display_genes == 'decreasing' and (which == 'euclid'):
                        idx = trg_sim.min(1).sort_values().index
                    elif which == 'intersect':
                        if display_genes.startswith('in') and mgt == 'down' or \
                        display_genes.startswith('de') and mgt == 'up':
                            asc = True
                        elif display_genes.startswith('in') and mgt == 'up' or \
                        display_genes.startswith('de') and mgt == 'down':
                            asc = False
                        idx = trg_sim.sum(1).sort_values(ascending=asc).index
                    get_genes = idx[:gene_number]
                    
                if specific_genes is not None:
                    # check if passed genelist in target markergenes add them 
                    # if not already in 
                    inp_ensg = util.get_ensgs(specific_genes, self._species)
                    not_mg = filter(lambda ie: ie not in trg_sim.index, inp_ensg)
                    inv = genes.set_index('ensg').loc[not_mg].name
                    if not inv.empty:
                        logger.info('{} not included: not markergenes of `'
                                    '{}-{}`'.format(inv.tolist(), mgt, trg))
                    add = lambda ie: not (ie in get_genes or ie in inv)
                    add_genes = pd.Index(filter(add, inp_ensg))
                    if not add_genes.empty:
                        get_genes = get_genes.append(add_genes)
                elif custom_target_genelist:
                    get_genes = genes.ensg
                
                if get_genes.empty:
                    logger.error('No genes were picked for {}-{}.'
                                 .format(mgt, trg))
                    sys.exit(1)
                # index overlap data to genelist
                ts = trg_sim.reindex(get_genes)
                ctrl_n = mgt, trg, ctrl_name

                if which == 'euclid':
                    cs = ctrl_sim.loc[get_genes, ctrl_n]
                    if not proportional:
                        agg_sim = ts.mean()
                    else:
                        s_np = s_noprop.loc[get_genes, (mgt, trg, slice(None))]
                        cs_np = cs_noprop.loc[get_genes, (mgt, trg, samples._ctrl)]
                        s_np = s_np.abs().mean()
                        cs_np = cs_np.abs().mean()
                        agg_sim = (cs_np - s_np) /cs_np
                elif which == 'intersect':
                    mgt_int = 1 if mgt == 'up' else -1
                    cs = pd.Series(mgt_int, get_genes, name=ctrl_n)
                    agg_sim = ts.sum() *mgt_int
                    if proportional:
                        agg_sim /= cs.shape
                # store plot data in nested dict
                data[trg][mgt] = (ts.T, cs.to_frame().T, agg_sim)
            
            # iterate target+markergene type
            sim.groupby(axis=1, level=(0,1), sort=False).apply(sel_genes, genes)
            return data

        # get data limits across all targets and markergene types to plot with 
        # one consistent heatmap range 
        def get_caps():
             # unpack nested dict into the 3 plot data elements       
            data_l = [e for dat in list(data.values())
                      for d in list(dat.values()) for e in d]
            tss, css, agss = [data_l[get::3] for get in (0,1,2)]
            
             # get number of genes per plot
            n_genes = [ts.shape[1] for ts in tss]
            if self._down_mgs:
                n_genes = [max(tss[i].shape[1], tss[i+1].shape[1]) 
                           for i in range(0, len(tss), 2)]

             # get sum plot limits
            agg_min = min([ags.min() for ags in agss])
            agg_max = max([ags.max() for ags in agss])
            agg_lim = [agg_min -abs(agg_min*.15), agg_max +abs(agg_max*.15)]
            # make sure 0 is included
            if differential:
                if agg_lim[0]>=0 and agg_lim[1]>=0:
                    agg_lim[agg_lim.index(min(agg_lim))] = 0
                elif agg_lim[0]<=0 and agg_lim[1]<=0:
                    agg_lim[agg_lim.index(max(agg_lim))] = 0
            if which == 'euclid':
                # get sum heatmap range
                mini = [ts.min().sort_values()[int(ts.shape[1]*.05)] for ts in tss]
                maxi = [ts.max().sort_values()[int(ts.shape[1]*.95)] for ts in tss]
                cap = round(max((abs(min(mini)), abs(max(maxi)))), 1)
                # get required_effect bar range
                re_cap = round(max([cs.iloc[0].sort_values()[int(cs.shape[1]*.95)]
                                   for cs in css]), 1)
                if not differential:
                    cap = re_cap = max((cap, re_cap))
                cap = 1 if proportional and (cap >1) else cap
                return cap, re_cap, agg_lim, n_genes
                
            if which == 'intersect' and not proportional:
                agg_lim = int(agg_lim[0]), int(agg_lim[1])
            return None, None, agg_lim, n_genes
        
        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            nplts = [4, 4]
            # default size of an exes is 0
            fig_widths = [.0001] *(nplts[1] +3)
            # based on parameters and config constants, set all sizes
            l = config.HM_LEFT if not pivot else config.HM_PIVOT_LEFT
            fig_widths[0] = samplelabels_space_left if samplelabels_space_left \
                            else l
            if show_samples_colorbar:
                fig_widths[1] = config.HM_Y_COLORBAR
            # heatmap width varies across plots, a nested list stores widths
            fig_widths[2] = [n_gs*config.HM_SQUARE_SIZE for n_gs in n_genes]
            if heatmap_width:
                fig_widths[2] = [heatmap_width*f_ws2 for f_ws2 in fig_widths[2]]
            if cluster_samples and show_samples_dendrogram:
                fig_widths[3] = config.HM_Y_DENDROGRAM
            if not hide_sum_plot:
                fig_widths[4] = config.G_HM_SUMPLOT_SIZE
            fig_widths[5] = config.HM_WSPACE * (nplts[1]-1)
            fig_widths[6] = samplelabels_space_right if samplelabels_space_right \
                            else config.HM_RIGHT

            fig_heights = [.0001] *(nplts[0] +3)
            fig_heights[0] = config.HM_TOP
            if cluster_genes and not hide_genes_dendrogram:
                fig_heights[1] = config.HM_X_DENDROGRAM
            if not hide_distance_bar:
                fig_heights[2] = config.HM_REQU_EFF_BAR 
            fig_heights[3] = config.HM_SQUARE_SIZE *len(samples._names_noctrl) 
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if show_genes_colorbar:
                fig_heights[4] = config.HM_X_COLORBAR
            fig_heights[5] = config.HM_HSPACE * (nplts[0]-1)
            b = config.HM_BOTTOM if not pivot else config.HM_PIVOT_BOTTOM
            fig_heights[6] = genelabels_space if genelabels_space else b

            # duplicate height sizes and insert a spacer axis with size of top
            if self._down_mgs:
                nplts[0] = nplts[0] *2 +1
                hs = fig_heights
                ins = [config.G_HM_UPDOWN_SPACE_SIZE]
                fig_heights = hs[:-2] + ins + hs[1:-2] + hs[-2:]
                fig_heights[-2] = config.HM_HSPACE * (nplts[0]-1)
            return nplts, fig_widths, fig_heights

        # draw plot
        def do_plot(i):
            # get final width list for specific number of genes in plot
            this_fig_widths = fig_widths[:2] +[fig_widths[2][i]] +fig_widths[3:]
            width, height = sum(this_fig_widths), sum(fig_heights)
            fig, axes = util._init_figure(this_fig_widths, fig_heights, nplts, 
                                          (config.HM_WSPACE, config.HM_HSPACE))
            if self._down_mgs:
                [ax.set_visible(False) for ax in axes[4, :]]

            # set plot title
            if title in ('None', 'none', 'False', 'false', 'F', 'f'):
                this_t = False
            if title:
                if title and type(title) is not str:
                    this_t = util._make_title(differential, proportional, which,
                                              samples.name, t_name, 
                                              postf='single gene ')
                else:
                    this_t = title
                if not pivot:
                    fig.suptitle(this_t, fontsize=config.FONTS *1.1,
                                 y=1-((config.HM_TOP/height *.15)))
                else:
                    axes[2, 0].set_ylabel(this_t, fontsize=config.FONTS *1.1,
                                            labelpad=18)
                
            # iterate over up and down plot-halfs
            for mgt, r in zip(self._mg_types, (0, 5)):
                sim, ctrl_sim, sim_agg = dat[mgt]

                # cluster genes/ samples and draw dendrograms
                if cluster_genes:
                    at = axes[r, 1] if not hide_genes_dendrogram else axes[r, 0]
                    order = util._heatmap_cluster(sim, 'top', at, 'columns')
                    sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
                if cluster_samples:
                    at = axes[2+r, 2] if show_samples_dendrogram else axes[r, 0]
                    order = util._heatmap_cluster(sim, 'right', at, 'rows')
                    sim, sim_agg = util._align_indices([sim, sim_agg], order, 0)
                axes[r, 0].set_visible(False)

                # draw the required_effect bar
                if not hide_distance_bar:
                    if reorder_to_distance_bar:
                        order = ctrl_sim.iloc[0].sort_values().index
                        sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
                    draw_cb = False
                    if which == 'euclid':
                        if not hide_colorbar_legend and (mgt == 'up') and \
                        differential:
                            draw_cb = True
                        cb_lbl = ('base expr. similarity\n'
                                  '[absolute eucl. distance]')
                        ctrl_lbl = samples._ctrl if not hide_samplelabels else ''
                        bar_args = {'vmin': 0, 'vmax': re_cap,
                                    'cmap': 'afmhot'}
                    elif which == 'intersect':
                        draw_cb = False
                        cb_lbl = None
                        ctrl_lbl = 'target markergenes' if not hide_samplelabels \
                                   else ''
                        bar_args = {'vmin': -1, 'vmax': 1,
                                    'cmap': config.RdBu_bin}
                    if distance_bar_range is not None:
                        bar_args.update({'vmin': distance_bar_range[0], 
                                         'vmax': distance_bar_range[1]})
                    util._plot_distance_bar(axes[1+r, :2], ctrl_sim, 
                                                  ctrl_lbl, bar_args, draw_cb, 
                                                  cb_lbl, fig, pivot, width, height)

                # setup heatmap x axis, including the colorbar
                xlbl = genes.set_index('ensg').reindex(sim.columns).name.values
                if show_genes_colorbar:
                    default = show_genes_colorbar.get('default', 'w') 
                    cols = [show_genes_colorbar.get(g, default) for g in xlbl]
                    cols = [c if is_color_like(c) else default for c in cols]
                else:
                    cols = None
                util._setup_heatmap_xy('x', axes[3+r, 1], xlbl, pivot,
                                      hide_genelabels, genelabels_size, cols) 


                # setup heatmap y axis, including the colorbar
                ylbl = sim.index.unique(2)[::-1]
                cols = samples.get_colors(ylbl) if show_samples_colorbar else \
                       None
                util._setup_heatmap_xy('y', axes[2+r, 0], ylbl, pivot, 
                                      hide_samplelabels, genelabels_size, cols)
                if self._down_mgs:
                    tit = '{} markergenes'.format(mgt)
                    axes[2+r, 0].set_title(tit, loc='right', fontweight='bold', 
                                           fontsize=config.FONTS)
                
                # draw sum plot on the right
                if not hide_sum_plot:
                    # general setup
                    ax = axes[2+r, 3]
                    ax.tick_params(labelbottom=True, bottom=True)
                    if pivot:
                        ax.tick_params(labelrotation=90)
                         
                    axes[3+r, 3].set_visible(False)
                    axes[1+r, 3].set_visible(False)
                    ax.set_axisbelow(True)
                    ax.xaxis.grid(alpha=0.8, linestyle='dashed')

                    # setup y axes
                    nsmps = sim_agg.shape[0]
                    ax.set_ylim(-.1, nsmps+.1)
                    yts = np.arange(nsmps-.5, -.5, -1)
                    ax.set_yticks(yts)
                    if show_samplelabels_sum_plot:
                        ax.tick_params(labelright=True)
                        ax.set_yticklabels(sim_agg.index.unique(2))

                    # setup x axes
                    xlim = agg_lim if not sum_plot_range else sum_plot_range
                    if pivot:
                        xlim.reverse()
                    ax.set_xlim(xlim)
                    if which == 'euclid' and differential and not proportional:
                        lbl = config.AGG_EUCL_DIFF_NOPROP
                    if which == 'euclid' and differential and proportional:
                        lbl = config.AGG_EUCL_DIFF_PROP
                    elif which == 'euclid' and not differential:
                        lbl = config.AGG_EUCL_NODIFF
                    elif which == 'intersect' and not proportional:
                        lbl = config.AGG_INTE_DIFF_NOPROP
                    elif which == 'intersect' and proportional:
                        lbl = config.AGG_INTE_DIFF_PROP
                    if (which == 'euclid') and not differential and samples._ctrl:
                        base = ctrl_sim.mean(1)
                        ax.vlines(base, 0, nsmps)
                        lbl = lbl[:26] + ' (line = base)' +lbl[26:]
                    ax.set_xlabel(lbl)
                        
                    if (which == 'euclid') and differential:
                        blue = config.colors[18] 
                        red = config.colors[14]
                        cols = [blue if v <0 else red for v in sim_agg.values]
                    else:
                        cols = config.colors[19]
                    ax.barh(y=yts, width=sim_agg, color=cols)
                
                # draw heatmap
                ax = axes[2+r, 1]
                ax.set_yticks(np.arange(0, sim.shape[0]))
                ax.set_xticks(np.arange(0, sim.shape[1]))
                if which == 'euclid':
                    if differential and not proportional:
                        hm_args = {'cmap': 'RdBu_r', 'vmin': -cap, 'vmax': cap}
                        cb_lbl = ('change in expr. similarity\n'
                                  '[differential eucl. distance]')
                    if differential and proportional:
                        hm_args = {'cmap': 'RdBu_r', 'vmin': -cap, 'vmax': cap}
                        cb_lbl = ('prop. of changed expr. similarity\n'
                                  '[prop. differential eucl. dist.]')
                    if not differential:
                        hm_args = {'cmap': 'afmhot', 'vmin': -0, 'vmax': cap}
                        cb_lbl = ('absolute expr. similarity\n'
                                  '[absolute eucl. distance]')
                elif which == 'intersect':
                    hm_args = {'cmap': config.RdBu_bin, 'vmin': -1, 'vmax': 1}
                    cb_lbl = ('differential genes similarity\n'
                              '[target markergene intersect]')
                
                if heatmap_range is not None:
                    hm_args.update({'vmin': heatmap_range[0], 
                                    'vmax': heatmap_range[1]})
                im = ax.imshow(sim.values, aspect='auto', **hm_args)

                # setup heatmap colorbar legend and draw
                if mgt == 'up' and not hide_colorbar_legend:
                    # add a new axis for the colorbar
                    at = (config.CB_LEFT/width, 1-config.CB_TOP/height, 
                          config.CB_WIDTH/width, config.CB_HEIGHT/height)
                    cax = fig.add_axes(at)
                    cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
    
                    cb.ax.set_xlabel(cb_lbl)
                    cb.ax.get_xaxis().set_label_position('top')
                    bar_ticks = [hm_args['vmin'], hm_args['vmax']]
                    cb.set_ticks(bar_ticks)
                    if which == 'intersect':
                        bar_ticks = ('down-regul.', 'up-regul.')
                    elif which == 'euclid' and proportional:
                        bar_ticks[0] = '< ' + str(bar_ticks[0])
                    cb.ax.set_xticklabels(bar_ticks)
                    if pivot:
                        cb.ax.tick_params(labelrotation=90)
                dat[mgt] = sim, ctrl_sim, sim_agg
            return fig, axes, dat
        
        spacer.info('\n\n' + log_plot)
        logger.info('Plot: {} & {}'.format(self.name, samples.name))
        genes = _check_args()
        data = get_data()
        cap, re_cap, agg_lim, n_genes = get_caps()

        nplts, fig_widths, fig_heights = get_plot_sizes()
        spacer.info('')
        logger.info('Drawing...')
        filename, pp = util._open_file(filename)
        ret = {}
        l_fn = filename.rfind('.png')
        for i, (t_name, dat) in enumerate(data.items()):
            fig, axes, dat = do_plot(i)
            spacer.info('{}/{} --- {}'.format(i+1, len(data), t_name))
            # plt.show()
            ret.update({t_name: (fig, axes, dat)})
            this_png_fn = '{}_{}.png'.format(filename[:l_fn], t_name)
            util._save_file(fig, filename=this_png_fn, pp=pp)
        if pp:
            pp.close()
        logger.info('Plots saved at {}/{}\n\n'
                    .format(os.path.abspath(os.curdir), filename))
        return ret 

    def ranked_similarity_barplot(self,
                                  # plot data
                                  samples,
                                  which = None, 
                                  differential = False,
                                  proportional = False,
                                  display_similarity = 'mgs mean',
                                  n_targets = 16,
                                  display_negative = True,
                                  # data ordering
                                  rank_samples = False,
                                  # general settings
                                  pivot = False,
                                  xlim_range = None,
                                  targetlabels_space = None,
                                  colored_bars = False,
                                  title = True,
                                  # show/ hide elements
                                  hide_targetlabels = False,
                                  hide_colorbar = False,
                                  filename = 'ranked_similarity_bp.pdf'):

        # check user input for errors and incompatibilities around `which` arg
        def _check_args():
            nonlocal which
            nonlocal differential
            nonlocal proportional
            nonlocal n_targets
            nonlocal display_similarity

            # check general basic input requirements
            a = util._check_args(self, samples, which, differential, proportional, 
                                display_similarity=display_similarity)
            which, differential, proportional, _, _, _, display_similarity = a
            if not n_targets:
                n_targets = len(self)
                logger.warning('The number of targets `n_targets` was not '
                               'passed. Set to all target elements ({}).'
                               .format(len(self)))
            logger.info('Arguments passed. Getting data now ...')

        # get the aggrevated overlap data for plotting, pick the targets
        def get_data():
            # dr_ctrl = True if differential else False
            _, sim, _, ctrl_sim = self._get_from_overlap(samples, which, 
                                                    genes_agg_diff=differential, 
                                                    genes_agg_prop=proportional,
                                                    )
            sim = sim.xs(display_similarity, 1, 0)
            if rank_samples:
                if differential:
                    order = sim.max(1).sort_values(ascending=False).index
                else:
                    order = sim.min(1).sort_values().index
                sim = sim.reindex(order)
            
            drop = slice(int(n_targets/2), -int(n_targets/2)) if display_negative \
                   else slice(-1, n_targets-1, -1)
            ascend = False if differential else True
            data = dict.fromkeys(sim.index, None)
            def sel_trgs(smp_row):
                trgs = smp_row.iloc[0].sort_values(ascending=ascend)
                data[trgs.name] = trgs.drop(trgs.index[drop])
            sim.groupby(level=0).apply(sel_trgs)
            return data, ctrl_sim
        
        # get plot global limits
        def get_caps():
            maxi = max([trg_vals.max() for trg_vals in list(data.values())])
            mini = min([trg_vals.min() for trg_vals in list(data.values())])
            lims = [mini -abs(mini*.15), maxi +abs(maxi*.15)]
            if lims[0]>=0 and lims[1]>=0:
                lims[lims.index(min(lims))] = 0
            elif lims[0]<=0 and lims[1]<=0:
                lims[lims.index(max(lims))] = 0
            return lims

        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            fig_widths = [.0001] *5
            l = config.BP_LEFT if not pivot else config.BP_PIVOT_LEFT
            fig_widths[0] = targetlabels_space if targetlabels_space else l
            if not hide_colorbar:
                fig_widths[1] = config.BP_Y_COLORBAR
            fig_widths[2] = config.BP_BARSPACE
            fig_widths[3] = .04
            r = config.BP_RIGHT if not pivot else config.BP_PIVOT_RIGHT
            fig_widths[4] = r

            fig_heights = [.0001] *4
            fig_heights[0] = config.BP_TOP
            fig_heights[1] = config.BP_BARWIDTH_SIZE *n_targets
            fig_heights[2] = 0
            fig_heights[3] = config.BP_BOTTOM
            return fig_widths, fig_heights
        
        # draw plot
        def do_plot(i, dat):
            height, width = sum(fig_heights), sum(fig_widths)
            fig, axes = util._init_figure(fig_widths, fig_heights, (1, 2), 
                                          (.04,0))
            ax = axes[1]

            # set plot title
            if title in ('None', 'none', 'False', 'false', 'F', 'f'):
                this_t = False
            if title:
                if title and type(title) is not str:
                    this_t = util._make_title(differential, proportional, which,
                                              s_name, self.name, pref='ranked ')
                else:
                    this_t = title
                if not pivot:
                    fig.suptitle(this_t, fontsize=config.FONTS *1.1, 
                                 y=1-(config.BP_TOP /height/6))
                else:
                    ax.get_yaxis().set_label_position('right')
                    ax.set_ylabel(this_t, rotation=-90, 
                                  fontsize=config.FONTS *1.1, labelpad=10)

            # setup y axis including the colorbar
            ax.spines['left'].set_visible(True)
            n = dat.shape[0] if not display_negative else dat.shape[0] +1
            ylim = n, -1
            yts = np.arange(n)
            [(ax.set_ylim(ylim), ax.set_yticks(yts)) for ax in axes]
            ylbls = dat.index.tolist()
            if not hide_colorbar:
                cols = self.get_colors(ylbls)
                if display_negative:
                    cols.insert(int(len(ylbls)/2), 'w')
                axes[0].bar(0, 1, color=cols, bottom=yts-.5)
            # if negative, insert a gab between the two groups 
            if display_negative:
                ylbls.insert(int(len(ylbls)/2), '')
                dat = dat.append(pd.Series(0, [''])).reindex(ylbls)
                # delta half-height/ width of split line between pos. & neg. group
                d_hh = (.01/fig_heights[1]) /2
                d_wh = (.03/fig_widths[2])
                line_args = {'xdata': (-d_wh, d_wh), 'transform': ax.transAxes, 
                             'clip_on': False, 'color': 'k'}
                ax.add_line(Line2D(ydata=(.5-d_hh*1.25, .5-d_hh*.25), **line_args))
                ax.add_line(Line2D(ydata=(.5+d_hh*.25, .5+d_hh*1.25), **line_args))
            if not hide_targetlabels:
                axes[0].tick_params(labelleft=True)
                if not pivot:
                    axes[0].set_yticklabels(ylbls)
                else:
                    axes[0].set_yticklabels(ylbls, rotation=-45, ha='right', 
                                            x=-.5, rotation_mode='anchor')
            
            # setup x axis
            xlim = lims
            if xlim_range:
                 xlim = xlim_range
            if not pivot:
                ax.spines['bottom'].set_visible(True)
                ax.tick_params(bottom=True, labelbottom=True)
            else:
                ax.spines['top'].set_visible(True)
                ax.tick_params(top=True, labeltop=True, labelrotation=-90)
                ax.xaxis.set_label_position('top')
            ax.set_xlim(xlim)
            ax.set_axisbelow(True)
            ax.xaxis.grid(alpha=0.8, linestyle='dashed')  
            if which == 'euclid' and differential and not proportional:
                xlbl = config.AGG_EUCL_DIFF_NOPROP
            if which == 'euclid' and differential and proportional:
                xlbl = config.AGG_EUCL_DIFF_PROP
            elif which == 'euclid' and not differential:
                xlbl = config.AGG_EUCL_NODIFF
            elif which == 'intersect' and not proportional:
                xlbl = config.AGG_INTE_DIFF_NOPROP
            elif which == 'intersect' and proportional:
                xlbl = config.AGG_INTE_DIFF_PROP
            # for absolute euclid sim., mark the untreated base if available
            if which == 'euclid' and not differential and samples._ctrl:
                xs = ctrl_sim.loc[samples._ctrl, display_similarity]
                xs = xs.reindex(ylbls, axis=1)
                ax.vlines(xs, yts-.4, yts+.4, linewidth=.5)
                xlbl = xlbl[:26] + ' (line = base)' +xlbl[26:]
            ax.set_xlabel(xlbl)

            if not colored_bars:
                cols = config.colors[19]
            else:
                blue = config.colors[18] 
                red = config.colors[14]
                cols = [blue if v <0 else red for v in dat.values]
            ax.barh(yts, dat, color=cols)
            return fig, axes

        spacer.info('\n\n' + log_plot)
        logger.info('Plot: {} & {}'.format(self.name, samples.name))
        _check_args()
        data, ctrl_sim = get_data()
        lims = get_caps()

        fig_widths, fig_heights = get_plot_sizes()
        spacer.info('')
        logger.info('Drawing...')
        filename, pp = util._open_file(filename)
        ret = {}
        l_fn = filename.rfind('.png')
        for i, (s_name, dat) in enumerate(data.items()):
            fig, axes = do_plot(i, dat)
            spacer.info('{}/{} --- {}'.format(i+1, len(data), s_name))
            # plt.show()
            ret.update({s_name: (fig, axes, dat)})
            this_png_fn = '{}_{}.png'.format(filename[:l_fn], s_name)
            util._save_file(fig, filename=this_png_fn, pp=pp)
        if pp:
            pp.close()
        logger.info('Plots saved at {}/{}\n\n'
                    .format(os.path.abspath(os.curdir), filename))
        return ret 