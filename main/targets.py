import pandas as pd
import numpy as np
import sys
import os
import copy

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import is_color_like


from DPre.main._differential import _differential
from DPre.main.drivers import Drivers
from DPre.main._logger import spacer, logger

import DPre.main.config as config
import DPre.main._dpre_util as util

class Targets(_differential):
    n_insts = 0
    def __init__(self, diff_genes, expression=None, ctrl=None, use_down_mgs=True, 
                 override_diff_names=False, name=None, log=True):
        
        self._down_mgs = use_down_mgs
        self._mg_types = ['up', 'down'] if use_down_mgs else ['up']
        super().__init__(diff_genes=diff_genes, expression=expression, ctrl=ctrl, 
                         override_diff_names=override_diff_names, name=name,
                         log=log, use_down_mgs=use_down_mgs)
        
        if self._ctrl:
            self._diff.drop(self._ctrl, axis=1, level=1, inplace=True)
            if self._has_expr:
                self._expr.drop(self._ctrl, axis=1, level=0, inplace=True)
            spacer.info('')
            logger.info('The control was dropped from `{}` since it has no '
                        'purpose in a Target instance.'.format(self._name))

        _diff_mgs = self._diff[self._diff.any(1)]
        self._diff_mgs = util._diff_to_int_updown_notation(_diff_mgs, False)

        if self._has_expr:
            expr_mgs = util._add_mg_types(self._expr.copy(), self._down_mgs)
            diff_masker = lambda trg: trg.mask(~self._diff[trg.columns[0][:-1]])
            expr_mgs = expr_mgs.groupby(level=(0,1), axis=1).apply(diff_masker)
            self._expr_mgs = expr_mgs.reindex(self._diff_mgs.index)

        self._overlaps = {}
        self.__class__.n_insts += 1


    def _overlap_drivers(self, drivers, which):
        logger.info('Overlappping Drivers ({}): `{}` on Targets: `{}` ... '
                    .format(which, drivers._name, self._name))

        if (which == 'euclid') and self._has_expr and drivers._has_expr:
            trg_data = self._expr_mgs.xs('z', 1, 2)
            drv_data = drivers._expr.xs('z', 1, 1)
            drv_data = util._add_mg_types(drv_data, self._down_mgs)
            subst = drivers._min_zval
        elif (which == 'intersect') and self._has_diff and drivers._has_diff:
            drv_data = drivers._diff_eff[self._mg_types]
            trg_data = self._diff_mgs
            trg_data.mask(trg_data == 0, inplace=True)
            subst = 0
        drv_data = drv_data.reindex(trg_data.index)

        drv_exts = [d for _, dd in drv_data.iteritems() 
                    for d in [dd]*len(self._names_noctrl)]
        drv_ext = pd.concat(drv_exts, axis=1)
        lvls = (self._mg_types, drv_ext.columns.unique(1), self._names_noctrl)
        drv_ext.columns = pd.MultiIndex.from_product(lvls).swaplevel()
        if config.UNDETECTED_MARKERGENES_BEHAVIOR == 'substitute':
            drv_ext.fillna(subst, inplace=True)
        elif config.UNDETECTED_MARKERGENES_BEHAVIOR == 'drop':
            drv_ext.dropna(inplace=True)
            trg_data = trg_data.reindex(drv_ext.index)

        eucl_dist = lambda drv_d: drv_d -trg_data[drv_d.columns.unique(0)].values 
        def mg_inters(drv_d):
            m = abs(drv_d + trg_data[drv_d.columns[0][0]].values)
            return m.mask(m == 2, 1).mask(m == 1, 0).mask(m == 0, -1)
        do_ovp = eucl_dist if which == 'euclid' else mg_inters
        ovp = drv_ext.groupby(axis=1, level=(0,2), sort=False).apply(do_ovp)
        

        self._overlaps['{}-{}'.format(id(drivers), which)] = ovp
        return ovp


    def _get_from_overlap(self, drivers, which, genes_diff=False, genes_prop=False,
                          norm_prop_genes=True,
                          genes_agg_diff=False, genes_agg_prop=False, 
                          idx_slice=pd.IndexSlice[:, :, :], 
                          drop_ctrl=True, inters_to_updown_not=True):

        # ==============================CHECK ARGS==============================
        if which == 'euclid':
            if idx_slice[2] != slice(None) and drivers._ctrl not in idx_slice[2]:
                idx_slice = (*idx_slice[:2], idx_slice[2]+[drivers._ctrl])
            if genes_agg_prop and not genes_agg_diff:
                genes_agg_prop = False
        try:
            ovp = self._overlaps['{}-{}'.format(id(drivers), which)]
        except KeyError:
            ovp = self._overlap_drivers(drivers, which)
        ovp = ovp.loc(1)[idx_slice].copy()

        if which == 'euclid':
            ctrl_genes = ovp.xs(drivers._ctrl, 1, 2, False)
            ctrl_agg = util.add_updown_mean(ctrl_genes.abs().mean().unstack().T)

            if not genes_diff:
                genes = ovp
            else:
                def gene_diff(drv_ovp, alter=None):
                    ctrl = drv_ovp.xs(drivers._ctrl, 1, 2).iloc(1)[0]
                    to_diff = lambda drv:  ctrl.abs() - drv.abs()
                    diff = drv_ovp.apply(to_diff)
                    if genes_prop:
                        to_prop = lambda drv: drv / ctrl.abs()
                        prop_interpr = ((diff.apply(to_prop)-1).abs() *-1) +1
                        return prop_interpr.mask(prop_interpr < -3, -3)
                    return diff
                genes = ovp.groupby(level=(0,1), axis=1, sort=False).apply(gene_diff)
               
            if not genes_agg_diff:
                agg = ovp.abs().mean().unstack((0,1)).reindex(ovp.columns.unique(2))
                agg = util.add_updown_mean(agg)
            else:
                def agg_diff(drv_ovp):
                    agg_mean = drv_ovp.abs().mean()
                    ctrl_mean = agg_mean.xs(drivers._ctrl, level=2).values
                    if genes_agg_diff:
                        agg_mean = ctrl_mean - agg_mean
                        if genes_agg_prop:
                            agg_mean /= ctrl_mean
                    return agg_mean.unstack((0,1)).reindex(agg_mean.index.unique(2))
                agg = ovp.groupby(level=(0,1), axis=1, sort=False).apply(agg_diff)
                agg = util.add_updown_mean(agg.droplevel((0,1), axis=1))
            
            if drop_ctrl:
                genes.drop(drivers._ctrl, axis=1, level=2, inplace=True)
                agg.drop(drivers._ctrl, inplace=True)
            return genes, agg, ctrl_genes, ctrl_agg
        
        elif which == 'intersect':
            if inters_to_updown_not and self._down_mgs:
                ovp['down'] *= -1
            n_mgs = ovp.notna().sum().xs(ovp.columns.unique(2)[0], level=2)
            n_mgs.index = util.add_level(n_mgs.index, 'n_mgs', at=2)
            n_mgs = util.add_updown_mean(n_mgs.unstack((0,1))).astype(int)

            agg = ovp.sum().unstack((0,1)).reindex(ovp.columns.unique(2))
            agg = util.add_updown_mean(agg)
            if genes_agg_prop:
                agg /= n_mgs.iloc[0]

            if drop_ctrl:
                ovp.drop(drivers._ctrl, axis=1, level=2, inplace=True)
                agg.drop(drivers._ctrl, inplace=True)
            return ovp, agg, None, n_mgs


    def target_similarity_heatmap(self, 
                                  drivers, 
                                  which, 
                                  proportional = False, 
   
                                  heatmap_range = None,
                                  heatmap_width = None,
                                  heatmap_height = None,
   
                                  show_required_effect_bar = True,
                                  reorder_to_required_effect_bar = False,
                                  required_effect_bar_range = None,
                                  
                                  cluster_targets = True,
                                  show_target_dendrogram = True,
                                  show_targets_colorbar = False,
                                  cluster_drivers = True,
                                  show_driver_dendrogram = True,
                                  show_drivers_colorbar = False,
                                  
                                  show_targetlabels = True,
                                  targetlabels_space = None,
                                  targetlabel_size = None,
                                  show_driverlabels = True,
                                  driverlabels_space = None,
                                  title = True, 
                                  show_colorbar_legend = True,
                                  filename = 'target_similarity_hm.png'):

        # check user input for errors and incompatibilities
        def check_args():
            nonlocal cluster_targets
            nonlocal reorder_to_required_effect_bar
            nonlocal show_required_effect_bar

            util.check_which(which, self, drivers)
            if reorder_to_required_effect_bar and not show_required_effect_bar:
                reorder_to_required_effect_bar = False
                logger.warning('When `reorder_to_required_effect_bar` is True, '
                               '`show_required_effect_bar` cannot be False. Set '
                               'to False.')
            if reorder_to_required_effect_bar and cluster_genes:
                cluster_genes = False
                logger.warning('Both `reorder_to_required_effect_bar` and '
                               '`cluster_genes` were set as True. '
                               '`cluster_genes` will be ignored.')

        # get the specific overlap data, plot the mean of up and down mgs
        def get_data():
            _, sim, _, ctrl_sim = self._get_from_overlap(drivers, which, 
                                                    genes_agg_diff=True, 
                                                    genes_agg_prop=proportional)
            return [sim.xs('mean', 1, 0), ctrl_sim.xs('mean', 1, 0)]

        # get plot lims
        def get_caps():
            mini = abs(data[0].min().min())
            maxi = abs(data[0].max().max())
            cap = round(max((mini, maxi)), 1)
            re_cap = round(data[1].iloc[0].max(), 1)
            return cap, re_cap
        
        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            nplts = [4,3]

            fig_widths = [.0001] *(nplts[1] +3)
            fig_widths[0] = driverlabels_space if driverlabels_space \
                            else config.TSH_LEFT
            if show_drivers_colorbar:
                fig_widths[1] = config.TSH_DRIVERS_COLORBAR
            fig_widths[2] = config.TSH_SQUARESIZE * len(self)
            if heatmap_width:
                fig_widths[2] *= heatmap_width
            if cluster_drivers and show_driver_dendrogram:
                fig_widths[3] = config.TSH_DRIVER_DENDROGRAM
            fig_widths[4] = config.TSH_WSPACE * (nplts[1]-1)
            fig_widths[5] = config.TSH_RIGHT

            fig_heights = [.0001] *(nplts[0] +3)
            fig_heights[0] = config.TSH_TOP
            if cluster_targets and show_target_dendrogram:
                fig_heights[1] = config.TSH_TARGET_DENDROGRAM
            if show_required_effect_bar:
                fig_heights[2] = config.TSH_REQU_EFF_BAR
            fig_heights[3] = config.TSH_SQUARESIZE * len(drivers._names_noctrl)
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if show_targets_colorbar:
                fig_heights[4] = config.TSH_TARGETS_COLORBAR
            fig_heights[5] = config.TSH_HSPACE * (nplts[0]-1)
            fig_heights[6] = targetlabels_space if targetlabels_space \
                             else config.TSH_BOTTOM

            return nplts, fig_widths, fig_heights

        # draw plot
        def do_plot():
            width, height = sum(fig_widths), sum(fig_heights)
            fig, axes = util._init_figure(fig_widths, fig_heights, nplts, 
                                       (config.TSH_WSPACE, config.TSH_HSPACE))
            sim, ctrl_sim = data

            # set plot title
            if title:
                if title and type(title) is not str:
                    this_title = ('transcriptinoal similarity\nof {} and {} '
                                  .format(drivers._name, self._name))
                fig.suptitle(this_title, fontsize=config.FONTS *1.2,
                             y=1-((config.GSH_TOP/height *.1)))


            # cluster targets/ drivers and draw dendrograms
            if cluster_targets:
                at = axes[0, 1] if show_target_dendrogram else axes[0, 0]
                order = util._heatmap_cluster(sim, 'top', at, 'columns')
                sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)
            if cluster_drivers:
                at = axes[2, 2] if show_driver_dendrogram else axes[0, 0]
                order = util._heatmap_cluster(sim, 'right', at, 'rows')
                # sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order, 0)
            axes[0, 0].set_visible(False)
            
            # draw required effect bar
            if show_required_effect_bar:
                if reorder_to_required_effect_bar:
                    order = ctrl_sim.iloc[0].sort_values().index
                    sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)

                ctrl_lbl = drivers._ctrl if show_driverlabels else ''
                draw_cb = True if show_colorbar_legend else False
                bar_args = {'cmap': 'afmhot', 'vmin': 0, 'vmax': re_cap}
                cb_lbl = 'number of markergenes\n' if which == 'intersect' else \
                         'base expr. similarity\n[mean abs. eucl. distance]'
                if required_effect_bar_range is not None:
                    bar_args.update({'vmin': required_effect_bar_range[0], 
                                     'vmax': required_effect_bar_range[1]})
                util.plot_required_effect_bar(axes[1, :2], ctrl_sim,
                                            ctrl_lbl, bar_args, draw_cb, 
                                            cb_lbl, fig, width, height)

            # setup heatmap x,y axis, including the colorbars
            util.setup_heatmap_xy('x', axes[3, 1], sim.columns, show_targetlabels, 
                                  targetlabel_size, show_targets_colorbar, 
                                  self.get_colors(sim.columns))
            util.setup_heatmap_xy('y', axes[2, 0], sim.index[::-1], show_driverlabels, 
                                  targetlabel_size, show_drivers_colorbar, 
                                  drivers.get_colors(sim.index[::-1]))

            # draw heatmap
            ax = axes[2, 1]
            ax.set_yticks(np.arange(0, sim.shape[0]))
            ax.set_xticks(np.arange(0, sim.shape[1]))
            hm_args = {'cmap': 'RdBu_r', 'vmin': -cap, 'vmax': cap}
            if which == 'euclid':
                if not proportional:
                    cb_lbl = config.AGG_EUCL_DIFF_NOPROP
                elif proportional:
                    cb_lbl = config.AGG_EUCL_DIFF_PROP
            elif which == 'intersect':
                if not proportional:
                    cb_lbl = config.AGG_INTE_DIFF_NOPROP
                elif proportional:
                    cb_lbl = config.AGG_INTE_DIFF_PROP
            if heatmap_range is not None:
                hm_args.update({'vmin': heatmap_range[0], 
                                'vmax': heatmap_range[1]})
            im = ax.imshow(sim.values, aspect='auto', **hm_args)
            
            # setup heatmap colorbar legend and draw
            if show_colorbar_legend:
                at = (config.CB_LEFT/width, 1-config.CB_TOP/height, 
                    config.CB_WIDTH/width, config.CB_HEIGHT/height)
                cax = fig.add_axes(at)
                cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
                
                bar_ticks = [hm_args['vmin'], hm_args['vmax']]
                if which == 'intersect' and not proportional:
                    bar_ticks = [int(bar_ticks[0]), int(bar_ticks[1])]
                cb.set_ticks(bar_ticks)
                cb.ax.set_xticklabels(bar_ticks)
                cb.ax.set_xlabel(cb_lbl)
                cb.ax.get_xaxis().set_label_position('top')

            return fig, axes, (sim, ctrl_sim)

        spacer.info('\n\n' + config.log_plot)
        check_args()
        logger.info('Processing data...')
        data = get_data()
        cap, re_cap = get_caps()

        nplts, fig_widths, fig_heights = get_plot_sizes()
        logger.info('Drawing plot: {} & {}'.format(self._name, drivers._name))
        pp = util.open_file(filename)
        fig, axes, data = do_plot()
        # plt.show()
        if filename.endswith('.pdf'):
            util.save_file(fig, pp)
            pp.close()
        elif filename.endswith('.png'):
            util.save_file(fig, filename)
        logger.info('Plot saved at {}'.format(os.path.abspath(os.curdir)))
        return fig, axes, data
        

    def gene_similarity_heatmap(self, 
                                drivers,  
                                which,
                                display_genes = 'increasing',
                                gene_number = 60,
                                specific_genes = None,
                                custom_target_genelist = None,
                                differential = True,
                                proportional = False, 
                                
                                heatmap_range = None,
                                heatmap_width = None,
                                heatmap_height = None,
            
                                show_required_effect_bar = True,
                                reorder_to_required_effect_bar = False,
                                required_effect_bar_range = None,
            
                                show_sum_plot = True,
                                sum_plot_range = None,
                                sum_plot_central_metric = 'mean',
                                show_driverlabels_sum_plot = False,
                                
                                cluster_genes = True,
                                show_gene_dendrogram = True,
                                genes_colorbar = None,
                                cluster_drivers = True,
                                show_driver_dendrogram = False,
                                show_drivers_colorbar = False,
                                
                                show_genelabels = True,
                                genelabel_space = None,
                                genelabel_size = None,
                                bold_emphz_genes = None,
                                show_driverlabels = True,
                                driverlabels_space = None,
                                title = True, 
                                show_colorbar_legend = True,
                                filename = 'gene_similarity_hm.pdf'):
        
        # check user input for errors and incompatibilities
        def check_args():
            nonlocal display_genes
            nonlocal specific_genes
            nonlocal custom_target_genelist
            nonlocal differential
            nonlocal proportional

            nonlocal show_required_effect_bar
            nonlocal reorder_to_required_effect_bar
            nonlocal sum_plot_central_metric
            nonlocal cluster_genes
            nonlocal genes_colorbar
        
            # check main data input
            util.check_which(which, self, drivers)
            if display_genes:
                val = ['variant', 'greatest', 'increasing', 'decreasing']
                if not differential:
                     val = ['variant', 'distant', 'similar']
                if display_genes not in val:
                    logger.error('The passed value for `display_genes`: {} is'
                                 'invalid. Valid options are {}'
                                 .format(display_genes, val))
                    sys.exit(1)
                if custom_target_genelist is not None:
                    custom_target_genelist = None
                    logger.info('Both `display_genes` and '
                                '`custom_target_genelist` were passed. '
                                '`custom_target_genelist` will be ignored.')
            elif custom_target_genelist is None and specific_genes is None:
                logger.error('None of `display_genes`, `specific_genes` or '
                             '`custom_target_genelist` were passed')
                sys.exit(1)
            elif custom_target_genelist is not None and specific_genes is not None:
                custom_target_genelist = None
                msg = ('Both `specific_genes` and `custom_target_genelist` were'
                       ' passed. `custom_target_genelist` will be ignored.')
                logger.info(msg)

            # inform user of smaller input incompatibilities 
            if which == 'intersect' and not differential:
                differential = True
                logger.warning('For the `intersect` similarity metric, '
                               'differential cannot be False. Was set to True.')
                if custom_target_genelist is not None and show_required_effect_bar:
                    show_required_effect_bar = False
                    msg = ('For the `intersect` similarity metric and a '
                           '`custom_target_genelist`, the `required_effect_bar`'
                           ' cannot be drawn. Set to False.')
                    logger.warning(msg)
            if proportional and not differential:
                proportional = False
                logger.warning('`proportional` can only be used if '
                               '`differential` is True aswell. Set to False.')
            if proportional and sum_plot_central_metric == 'mean':
                sum_plot_central_metric = 'median'
                msg = ('`proportional` and `sum_plot_central_metric` = mean'
                       'is not meaningful because outlaying negative values '
                       'are capped at -3. The median is used instead.')
                logger.warning(msg)
            if reorder_to_required_effect_bar and not show_required_effect_bar:
                reorder_to_required_effect_bar = False
                logger.warning('When `reorder_to_required_effect_bar` is True, '
                               '`show_required_effect_bar` cannot be False. Set '
                               'to False.')

            if reorder_to_required_effect_bar and cluster_genes:
                cluster_genes = False
                logger.warning('Both `reorder_to_required_effect_bar` and '
                               '`cluster_genes` were set as True. '
                               '`cluster_genes` will be ignored.')
            
            # modify arguments for convneience 
            if which == 'euclid' and not differential:
                if display_genes == 'distant':
                    display_genes = 'increasing'
                elif display_genes == 'similar':
                    display_genes = 'decreasing'
            if genes_colorbar == True:
                if specific_genes:
                    genes_colorbar = dict.fromkeys(specific_genes, 
                                                   config.colors[1])
                else:
                    genes_colorbar = None

            # get a list of generally valid annotated genes
            genes = pd.DataFrame({'name': util.annotate(self._diff_mgs.index), 
                                  'ensg': self._diff_mgs.index })
            if specific_genes is not None or custom_target_genelist is not None:
                # for gene input check if genes are detected in the target data
                if specific_genes is not None:
                    inp_gl = pd.Index(specific_genes).drop_duplicates()
                    val_gl = pd.Index(genes.name.values)
                    isin = 'markergenes'
                elif custom_target_genelist is not None:
                    inp_gl = pd.Index(custom_target_genelist).drop_duplicates()
                    val_gl = self._detec_genes
                    isin = 'detected genes'
                    # if the overlap contains only driver-detected genes:
                    if config.UNDETECTED_MARKERGENES_BEHAVIOR == 'drop':
                        val_gl = val_gl.intersection(drivers._detec_genes)
                        isin = 'valid genes'
                    val_gl = pd.Index(util.annotate(val_gl))
                
                inv = [g for g in inp_gl if g not in val_gl]
                inp_gl = inp_gl.drop(inv) 
                if inv:
                    logger.warning('{} ({}/{}) are not {} in any of the targets'
                                   ' in `{}`. These genes will not be included.'
                                   .format(inv, len(inv), len(inv)+len(inp_gl), 
                                           isin, self._name))
                    if len(inv) == (len(inv)+len(inp_gl)):
                        sys.exit(1)
                # update passed list
                if specific_genes is not None:
                    specific_genes = inp_gl
                elif custom_target_genelist is not None:
                    genes = pd.DataFrame({'name': inp_gl, 
                                          'ensg': util.get_ensgs(inp_gl)})
            return genes                   

        # get the specific overlap data and pick out the genes to display
        def get_data():
            # init a new target where all genes are markergenes all targets
            if custom_target_genelist:
                nonlocal self
                diff = pd.DataFrame(True, genes.ensg, self._diff.columns)
                args = {'diff_genes': diff}
                if which == 'euclid':
                    expr = self._expr.loc[genes.ensg].copy()
                    args.update({'expression': expr})
                self = Targets(name='custom genelist', use_down_mgs=False, 
                               log=False, **args)
            sim, _, ctrl_sim, _ = self._get_from_overlap(drivers, which, 
                                                      genes_diff=differential, 
                                                      genes_prop=proportional, 
                                                      norm_prop_genes=True)
            if which == 'euclid':
                ctrl_sim = ctrl_sim.abs()
                if not differential:
                    sim = sim.abs()

            # init mutable dict with target and markegene type keys
            data = dict((trg, dict((mgt, None) for mgt in self._mg_types))
                        for trg in self._names_noctrl)
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
                    inp_ensg = util.get_ensgs(specific_genes)
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
                
                # index overlap data to genelist
                if get_genes.empty:
                    logger.error('No genes were picked for {}-{}.'
                                 .format(mgt, trg))
                    sys.exit(1)
                ts = trg_sim.reindex(get_genes)
                ctrl_name = mgt, trg, drivers._ctrl
                if which == 'euclid':
                    cs = ctrl_sim.loc[get_genes, ctrl_name]
                elif which == 'intersect':
                    mgt_int = 1 if mgt == 'up' else -1
                    cs = pd.Series(mgt_int, get_genes, name=ctrl_name)
                    agg_sim = ts.sum() *mgt_int
                    if proportional:
                        agg_sim /= cs.shape
                if (which == 'euclid'):
                    if sum_plot_central_metric == 'mean':
                        agg_sim = ts.mean()
                    elif sum_plot_central_metric == 'median':
                        agg_sim = ts.median()
                # store plot data in nested dict
                data[trg][mgt] = (ts.T, cs.to_frame().T, agg_sim)
            
            # iterate target+markergene type
            spacer.info('') 
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
            fig_widths[0] = driverlabels_space if driverlabels_space \
                            else config.GSH_LEFT
            if show_drivers_colorbar:
                fig_widths[1] = config.GSH_DRIVERS_COLORBAR
            # heatmap width varies across plots, a nested list stores widths
            fig_widths[2] = [n_gs*config.GSH_SQUARESIZE for n_gs in n_genes]
            if heatmap_width:
                fig_widths[2] = [heatmap_width*f_ws2 for f_ws2 in fig_widths[2]]
            if cluster_drivers and show_driver_dendrogram:
                fig_widths[3] = config.GSH_DRIVERS_DENDROGRAM
            if show_sum_plot:
                fig_widths[4] = config.GSH_SUMPLOT_SIZE
            fig_widths[5] = config.GSH_WSPACE * (nplts[1]-1)
            fig_widths[6] = config.GSH_RIGHT

            fig_heights = [.0001] *(nplts[0] +3)
            fig_heights[0] = config.GSH_TOP
            if cluster_genes and show_gene_dendrogram:
                fig_heights[1] = config.GSH_GENE_DENDROGRAM
            if show_required_effect_bar:
                fig_heights[2] = config.GSH_REQU_EFF_BAR
            fig_heights[3] = config.GSH_SQUARESIZE *len(drivers) 
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if genes_colorbar:
                fig_heights[4] = config.GSH_GENES_COLORBAR
            fig_heights[5] = config.GSH_HSPACE * (nplts[0]-1)
            fig_heights[6] = genelabel_space if genelabel_space \
                             else config.GSH_BOTTOM

            # duplicate height sizes and insert a spacer axis with size of top
            if self._down_mgs:
                nplts[0] = nplts[0] *2 +1
                hs = fig_heights
                fig_heights = hs[:-2]+[config.GSH_UP_DOWN_SPACE]+hs[1:-2]+hs[-2:]
                fig_heights[-2] = config.GSH_HSPACE * (nplts[0]-1)
            return nplts, fig_widths, fig_heights

        # draw plot
        def do_plot(i):
            # get final width list for specific number of genes in plot
            this_fig_widths = fig_widths[:2] +[fig_widths[2][i]] +fig_widths[3:]
            width, height = sum(this_fig_widths), sum(fig_heights)
            fig, axes = util._init_figure(this_fig_widths, fig_heights, nplts, 
                                          (config.GSH_WSPACE, config.GSH_HSPACE))
            if self._down_mgs:
                [ax.set_visible(False) for ax in axes[4, :]]

            # set plot title
            if title:
                if title and type(title) is not str:
                    this_title = ('single-gene transcriptional similarity\nof '
                                 '{} and {} '.format(drivers._name, t_name))
                else:
                    this_title = title
                fig.suptitle(this_title, fontsize=config.FONTS *1.2,
                             y=1-((config.GSH_TOP/height *.15)))
                
            # iterate over up and down plot-halfs
            for mgt, r in zip(self._mg_types, (0, 5)):
                sim, ctrl_sim, sim_agg = dat[mgt]

                # cluster genes/ drivers and draw dendrograms
                if cluster_genes:
                    at = axes[r, 1] if show_gene_dendrogram else axes[r, 0]
                    order = util._heatmap_cluster(sim, 'top', at, 'columns')
                    sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)
                if cluster_drivers:
                    at = axes[2+r, 2] if show_driver_dendrogram else axes[r, 0]
                    order = util._heatmap_cluster(sim, 'right', at, 'rows')
                    sim, sim_agg = util.align_indices([sim, sim_agg], order, 0)
                axes[r, 0].set_visible(False)

                # draw the required_effect bar
                if show_required_effect_bar:
                    if reorder_to_required_effect_bar:
                        order = ctrl_sim.iloc[0].sort_values().index
                        sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)
                    draw_cb = False
                    if which == 'euclid':
                        if show_colorbar_legend and (mgt == 'up') and differential:
                            draw_cb = True
                        cb_lbl = ('base expr. similarity\n'
                                  '[absolute eucl. distance]')
                        ctrl_lbl = drivers._ctrl if show_driverlabels else ''
                        bar_args = {'vmin': 0, 'vmax': re_cap,
                                    'cmap': 'afmhot'}
                    elif which == 'intersect':
                        draw_cb = False
                        cb_lbl = None
                        ctrl_lbl = 'target markergenes' if show_driverlabels else ''
                        bar_args = {'vmin': -1, 'vmax': 1,
                                    'cmap': config.RdBu_bin}
                    if required_effect_bar_range is not None:
                        bar_args.update({'vmin': required_effect_bar_range[0], 
                                         'vmax': required_effect_bar_range[1]})
                    util.plot_required_effect_bar(axes[1+r, :2], ctrl_sim, 
                                                  ctrl_lbl, bar_args, draw_cb, 
                                                  cb_lbl, fig, width, height)

                # setup heatmap x axis, including the colorbar
                xlbl = genes.set_index('ensg').reindex(sim.columns).name.values
                if genes_colorbar:
                    default = genes_colorbar.get('default', 'w') 
                    cols = [genes_colorbar.get(g, default) for g in xlbl]
                    cols = [c if is_color_like(c) else default for c in cols]
                else:
                    cols = None
                util.setup_heatmap_xy('x', axes[3+r, 1], xlbl, show_genelabels, 
                                      genelabel_size, genes_colorbar, cols) 

                # setup heatmap y axis, including the colorbar
                ylbl = sim.index.unique(2)[::-1]
                util.setup_heatmap_xy('y', axes[2+r, 0], ylbl, show_driverlabels, 
                                      genelabel_size, show_drivers_colorbar, 
                                      drivers.get_colors(ylbl))
                if self._down_mgs:
                    tit = '{} markergenes'.format(mgt)
                    axes[2+r, 0].set_title(tit, loc='right', fontweight='bold', 
                                           fontsize=config.FONTS)

                # draw sum plot on the right
                if show_sum_plot:
                    # general setup
                    metric = sum_plot_central_metric
                    ax = axes[2+r, 3]
                    ax.tick_params(labelbottom=True, bottom=True)
                    axes[3+r, 3].set_visible(False)
                    axes[1+r, 3].set_visible(False)
                    ax.set_axisbelow(True)
                    ax.xaxis.grid(alpha=0.8, linestyle='dashed')

                    # setup x axes
                    xlim = agg_lim if not sum_plot_range else sum_plot_range
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
                    if sum_plot_central_metric == 'median':
                        lbl = lbl.replace('mean', 'median')
                    ax.set_xlabel(lbl)
                    
                    # setup y axes
                    ndrvs = sim_agg.shape[0]
                    ax.set_ylim(-.1, ndrvs+.1)
                    yts = np.arange(ndrvs-.5, -.5, -1)
                    ax.set_yticks(yts)
                    if show_driverlabels_sum_plot:
                        ax.tick_params(labelright=True)
                        ax.set_yticklabels(sim_agg.index.unique(2))

                    if (which == 'euclid') and not differential:
                        base = ctrl_sim.mean(1) if metric == 'mean' \
                               else ctrl_sim.median(1)
                        ax.vlines(base, 0, ndrvs)
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
                if mgt == 'up' and show_colorbar_legend:
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
                dat[mgt] = sim, ctrl_sim, sim_agg
            return fig, axes, dat
      
      
        spacer.info('\n\n' + config.log_plot)
        genes = check_args()
        logger.info('Processing data...')
        data = get_data()
        spacer.info('') 
        cap, re_cap, agg_lim, n_genes = get_caps()

        logger.info('Drawing plot: {} & {}'.format(self._name, drivers._name))
        nplts, fig_widths, fig_heights = get_plot_sizes()
        filename, pp = util.open_file(filename)
        ret = {}
        for i, (t_name, dat) in enumerate(data.items()):
            fig, axes, dat = do_plot(i)
            logger.info('Drawing plot: {} & {}'.format(t_name, drivers._name))
            # plt.show()
            ret.update({t_name: (fig, axes, dat)})
            if filename.endswith('.pdf'):
                util.save_file(fig, filename=filename, pp=pp)
            elif filename.endswith('.png'):
                util.save_file(fig, filename=filename[:-4] +'_' +t_name +'.png')
        if filename.endswith('.pdf'):
            pp.close()
        logger.info('Plots saved at {}'.format(os.path.abspath(os.curdir)))
        return ret 
    
















    # DEPRECATED
    # DEPRECATED
    # DEPRECATED
    # DEPRECATED
    # DEPRECATED
    # DEPRECATED
    # DEPRECATED

    def ranked_similarity_barplot(self,
                                  drivers,
                                  which, 
                                  plot_celltype_identity = False,
                                  proportional = False,
                                  target_numbers = None,
                                  target_threshhold = None,
                                  driver_numbers = None,
      
                                  show_updown = False,
                                  plot_layout = None,
                                  xlim_from_config = False,
                                  load_overlap = True,
                                  labelsize = 1,
                                  show_target_labels = True,
                                  ylabel_space = None,
                                  colored = False,
                                  title = True,
                                  colorbar = True,
                                  filename = 'driver_effects2.pdf'):

        # ==============================CHECK ARGS==============================
        if show_updown and not self._down_mgs:
            show_updown = False
            logger.warning('No down markergenes in {}. `show_updown` == True '
                            'is invalid. Set to False.'.format(self.name))
        if not driver_numbers:
            driver_numbers = len(drivers)
        if not target_numbers:
            target_numbers = len(self)
        asc = True if plot_celltype_identity else False

        # =================================DATA=================================
        data = self.get_from_overlap(drivers, which, eff_prop=proportional, drop_ctrl=False,
                                     standardize_down_mgs=True, eff_diff=not asc)
        # get sigle effects, slice with passed arguaments to relevent ones
        data = data[1].reset_index(level=0, drop=True)
        if plot_celltype_identity:
            sor = data.xs('mean', level=1).min(axis=1).sort_values(ascending=asc)
        else:    
            sor = data.xs('mean', level=1).max(axis=1).sort_values(ascending=asc)
            
        data = data.reindex(sor.index[:driver_numbers], level=0)
        data[data < 0] = np.nan
        # print(data)
        data = data.groupby(level=0).apply(util.filter_trgs, target_numbers, 
                                           target_threshhold, ascending=asc)

        # =================================PLOT=================================
        pp = util.open_file(filename)
        drv_ns = data.index.unique(0)

        # set plotlayout if not passed to 3 columns
        if not plot_layout:
            if not len(drv_ns) %3:
                plot_layout = (int(len(drv_ns)/3), 3)
            else:
                plot_layout = (int(len(drv_ns)/3) +1, 3)
        nplts = plot_layout[0], plot_layout[1] *3
        
        # get max targets per plot and plot-row to adjust height accordingly
        n_trgs = data.xs('mean', level=1).notna().sum(axis=1)
        n_cols = plot_layout[1]
        row_n_trgs = [max(n_trgs[i*n_cols: i*n_cols +n_cols]) 
                      for i in range(len(n_trgs)) if i*n_cols < len(n_trgs)]
        
        # `bar_ws` for bar width in inches for defining figure height
        bar_ws = config.RSB_BARWIDTH_SIZE + config.RSB_BETWEEN_BARS_SIZE
        if show_updown:
            bar_ws += config.RSB_BARWIDTH_SIZE *.5
        # `btw_bars` and `bar_width` define sizes in plot scale
        btw_bars = config.RSB_BETWEEN_BARS_SIZE /config.RSB_BARWIDTH_SIZE
        bar_width = 1 +btw_bars if not show_updown else 1.5 +btw_bars

        # set widths in inches for every axis element
        fig_widths = [.0001] *5
        fig_widths[0] = config.RSB_LEFT
        fig_widths[1] = ylabel_space if ylabel_space else config.RSB_YLABELSPACE
        if colorbar:
            fig_widths[2] = config.RSB_TARGETS_COLORBAR
        fig_widths[3] = config.RSB_BARSPACE
        fig_widths[4] = config.RSB_RIGHT
        f_ws = fig_widths
        fig_widths = [f_ws[0]] + f_ws[1:-1]*plot_layout[1] + [f_ws[-1]]
        fig_widths.insert(-1, config.RSB_WSPACE * (nplts[1]-1))

        # set heights in inches for every axis element and the space between
        fig_heights = [.0001] *(plot_layout[0] +2)
        fig_heights[0] = config.RSB_TOP
        fig_heights[1:-1] = [bar_ws *n for n in row_n_trgs]
        fig_heights[-1] = config.RSB_BOTTOM
        fig_heights.insert(-1, config.RSB_HSPACE * (nplts[0]-1))

        # init barplot with set size
        width, height = sum(fig_widths), sum(fig_heights)
        fig, axes = util._init_figure(fig_widths, fig_heights, nplts, 
                                       (config.RSB_WSPACE, config.RSB_HSPACE))
        axes = axes.flatten()

        # main title
        if title:
            if title and type(title) is not str:
                if plot_celltype_identity:
                    title = ('{} cell type identity: {}'
                            .format(drivers.name, self.name))
                else:
                    title = ('{} target prediction: {}'
                            .format(drivers.name, self.name))

            fig.suptitle(title, fontsize=config.FONTS *1.2, 
                         y=1-(config.RSB_TOP /height/2))

        # iterate over the single plots
        row = 0
        for i in range(len(drv_ns)):
            axe = axes[i*3:i*3 +3]
            dat = data.loc[drv_ns[i]].dropna(axis=1)
            if i and not i % plot_layout[1]:
                row +=1

            # set basics for main plot like spines, grid and title
            ax = axe[2]
            ax.spines['bottom'].set_visible(True)
            ax.spines['left'].set_visible(True)
            ax.set_axisbelow(True)
            ax.xaxis.grid(alpha=0.8, linestyle='dashed')  
            ax.set_title(drv_ns[i], fontsize=config.FONTS *1.1, pad=1)

            # y-axis lim, ticks labels (linked to colorbar axes)
            c_yts = np.arange(0, dat.shape[1]*bar_width, bar_width)
            # occasional floating inaccuracy
            if len(c_yts) != dat.shape[1]:
                c_yts = c_yts[:dat.shape[1]]
            ylim = row_n_trgs[row]*bar_width, -bar_width/2
            axe[1].set_ylim(ylim)
            axe[2].set_ylim(ylim)
            axe[1].set_yticks(c_yts)
            axe[2].set_yticks(c_yts)
            ylbls = dat.sort_values(by='mean', ascending=asc, axis=1).columns
            if show_target_labels:
                fs = config.FONTS *labelsize
                axe[1].tick_params(left=True, labelleft=True)
                axe[1].set_yticklabels(ylbls, fontsize=fs)

            # x axis setup
            if not xlim_from_config:
                xlim = data.max().max() + data.max().max()/8
            else:
                 xlim = config.RSB_XLIM
            ax.tick_params(bottom=True, labelbottom=True)
            ax.set_xlim(0, xlim)
            # only show xlabel for plots in last respective row
            if i >= len(data.index.unique(0)) -plot_layout[1]:
                if which == 'euclid':
                    xlbl = 'z-trans. eucliden distance'
                    if not plot_celltype_identity:
                        xlbl += ' reduction'
                elif which == 'intersect':
                    xlbl = 'differential gene intersection'
                if proportional:
                    xlbl = 'proportional  ' + xlbl
                ax.set_xlabel(xlbl)

            # setup colorbar indicating targets
            if colorbar:
                bar_args = config.color_bar()
                bar_width_no_btw = bar_width -btw_bars
                tick_adj = 0 if not show_updown else .25
                bottom = [t +tick_adj -bar_width_no_btw/2 for t in c_yts]
                
                colors = self.get_colors(ylbls)
                for j in range(len(colors)):
                    axe[1].bar(0, bar_width_no_btw, color=colors[j], 
                               bottom=bottom[j], **bar_args)
            # draw barplot
            dat = dat.sort_values(by='mean', ascending=asc, axis=1)
            bar_args = config.default_single_hbarplot(0)
            if colored:
                cols = self.get_colors(ylbls)
                bar_args['color'] = cols
            axe[2].barh(y = c_yts, 
                        width = dat.loc['mean'], 
                        height = 1, 
                        **bar_args)

            # up- and down colorbar
            if show_updown:
                up = dat.loc['up'] /2
                down = dat.loc['down'] /2
                down_neg = (up >0) & (down <0)
                if any(down_neg.values):
                    up.loc[down_neg.values] = dat.loc['mean']
                    down.loc[down_neg.values] = 0
                
                up_args = config.default_single_hbarplot(1)
                ax.barh(y = c_yts +.75, 
                        width = up, 
                        height = .5, 
                        label = 'up-regul. genes', 
                        **up_args)
                down_args = config.default_single_hbarplot(2)
                ax.barh(y = c_yts +.75,
                        width = down, 
                        height = .5, 
                        left = up, 
                        label = 'down-regul. genes', 
                        **down_args)
                if i == 0:
                    # up and down bar legend
                    fig.legend(ncol=2, loc='lower center', 
                               bbox_to_anchor=(.5, config.RSB_BOTTOM/height *1/10))
        # plt.show()
        util.save_file(fig, filename, pp, last=True)
