import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches


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
        self._diff_mgs = util._diff_to_int_updown_notation(_diff_mgs)

        if self._has_expr:
            expr_mgs = util._add_mg_types(self._expr.copy(), self._down_mgs)
            diff_masker = lambda trg: trg.mask(~self._diff[trg.columns[0][:-1]])
            expr_mgs = expr_mgs.groupby(level=(0,1), axis=1).apply(diff_masker)
            self._expr_mgs = expr_mgs.reindex(self._diff_mgs.index)

        self._overlaps = {}
        self.__class__.n_insts += 1


    def _overlap_drivers(self, drivers, which):
        logger.info('Overlappping ({}) Drivers:`{}` on Targets: `{}` ... '
                    .format(which, drivers._name, self._name))

        if (which == 'euclid') and self._has_expr and drivers._has_expr:
            trg_data = self._expr_mgs.xs('z', 1, 2)
            drv_data = drivers._expr.xs('z', 1, 1)
            drv_data = util._add_mg_types(drv_data, self._down_mgs)
            subst = drivers._min_zval
            
        elif (which == 'intersect') and self._has_diff and drivers._has_diff:
            trg_data = util._add_mg_types(self._diff_mgs, self._down_mgs)
            trg_data.mask(trg_data == 0, inplace=True)
            drv_data = util._add_mg_types(drivers._diff_eff, self._down_mgs)
            subst = 0

        drv_data = drv_data.reindex(trg_data.index)
        drv_exts = [d for _, dd in drv_data.iteritems() for d in [dd]*len(self)]
        drv_ext = pd.concat(drv_exts, axis=1)
        lvls = (self._mg_types, drv_ext.columns.unique(1), self._names)
        drv_ext.columns = pd.MultiIndex.from_product(lvls).swaplevel()

        if config.UNDETECTED_MARKERGENES_BEHAVIOR == 'substitute':
            drv_ext.fillna(subst, inplace=True)
        elif config.UNDETECTED_MARKERGENES_BEHAVIOR == 'drop':
            drv_ext.dropna(inplace=True)
            trg_data = trg_data.reindex(drv_ext.index)

        eucl_dist = lambda drv_d: drv_d -trg_data[drv_d.columns.unique(0)].values 
        def mg_inters(drv_d):
            trg_d = trg_data[drv_d.columns[0][0]]
            return drv_d.mask(drv_d == trg_d.values, 1).mask(trg_d.isna().values)
        do_ovp = eucl_dist if which == 'euclid' else mg_inters
        ovp = drv_ext.groupby(axis=1, level=(0,2), sort=False).apply(do_ovp)

        self._overlaps['{}-{}'.format(id(drivers), which)] = ovp
        return ovp



    def _get_from_overlap(self, drivers, which, genes_diff=False, genes_prop=False,
                          norm_prop_genes=True,
                          genes_agg_diff=False, genes_agg_prop=False, 
                          idx_slice=pd.IndexSlice[:, :, :], 
                          drop_ctrl=True, get_max_possible_only=False):

        # ==============================CHECK ARGS==============================
        if which == 'euclid':
            if idx_slice[2] != slice(None) and drivers._ctrl not in idx_slice[2]:
                idx_slice = (*idx_slice[:2], idx_slice[2]+[drivers._ctrl])
            if get_max_possible_only and genes_agg_diff:
                genes_agg_diff = False
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
            if get_max_possible_only:
                return ctrl_genes, ctrl_agg

            if not genes_diff:
                genes = ovp
            else:
                def gene_diff(drv_ovp, alter=None):
                    ctrl = drv_ovp.xs(drivers._ctrl, 1, 2).iloc(1)[0]
                    to_diff = lambda drv:  ctrl.abs() - drv.abs()
                    diff = drv_ovp.apply(to_diff)
                    if genes_prop:
                        to_prop = lambda drv: drv / ctrl.abs()
                        diff = ((diff.apply(to_prop)-1).abs() *-1) +1
                        # if norm_prop_genes:
                        #     cpv = config.NORM_PROP_BAHVIOuR['change_prop_value']
                        #     alt = alter.xs(diff.columns[0][1], 1, 1, False)
                        #     alt_mask = np.repeat(alt.notna().values, diff.shape[1], 1)
                        #     if cpv == 'min_found':
                        #         subst = diff[alt.isna()].min().min()
                        #     elif isinstance(cpv, (int, float)):
                        #         subst = cpv
                        #     elif cpv == 'drop':
                        #         subst = np.nan
                        #     elif cpv == 'none':
                        #         return diff
                        #     diff = diff.mask((alt_mask), subst)
                        #     # plt.hist(diff.unstack().dropna(), density=True, bins=400, color='green', alpha=.6)
                        #     return diff
                    return diff

                genes = ovp.groupby(level=(0,1), axis=1, sort=False).apply(gene_diff)
                
                # genes_grp = ovp.groupby(level=(0,1), axis=1, sort=False)
                # if not (genes_prop and norm_prop_genes):
                #     genes = genes_grp.apply(gene_diff)
                # else:
                #     cap_type = config.NORM_PROP_BAHVIOuR['cap_type']
                #     cap = config.NORM_PROP_BAHVIOuR['cap']
                #     ctrls = ovp.xs(drivers._ctrl, 1, 2, False).abs()
                #     cumsum, vals = plt.hist(ctrls.dropna().unstack(), bins=250,
                #                             cumulative=1, density=True,)[:-1]
                #     plt.close()
                #     if cap_type == 'proportion':
                #         cap = vals[len(cumsum[cumsum <= cap])]
                #     alter = ctrls.mask(ctrls>cap)
                #     genes = genes_grp.apply(gene_diff, alter)

                    # print(len(genes.unstack().dropna().values))
                    # plt.hist(genes.unstack().dropna().values, density=False, bins=300)
                    # plt.show()


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

            n_mgs = ovp.notna().sum().xs(ovp.columns.unique(2)[0], level=2)
            n_mgs.index = util.add_level(n_mgs.index, 'n_mgs', at=2)
            n_mgs = util.add_updown_mean(n_mgs.unstack((0,1))).astype(int)
            if get_max_possible_only:
                return None, n_mgs

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
   
                                  show_max_possible_bar = True,
                                  reorder_to_max_possible_bar = False,
                                  max_possible_bar_range = None,
                                  
                                  cluster_targets = True,
                                  show_target_dendrogram = True,
                                  targets_colorbar = False,
                                  cluster_drivers = True,
                                  show_driver_dendrogram = True,
                                  drivers_colorbar = False,
                                  
                                  show_targetlabels = True,
                                  targetlabels_space = None,
                                  targetlabel_size = None,
                                  show_driverlabels = True,
                                  driverlabels_space = None,
                                  title = True, 
                                  show_colorbar_legend = True,
                                  filename = 'target_similarity_hm.png'):
          
        def check_args():
            if reorder_to_max_possible_bar:
                cluster_targets = False
            assert which in ('euclid', 'intersect'), '`which` invalid '
        
        def get_data():
            _, sim, _, ctrl_sim = self._get_from_overlap(drivers, which, 
                                                    genes_agg_diff=True, 
                                                    genes_agg_prop=proportional)
            return [sim.xs('mean', 1, 0), ctrl_sim.xs('mean', 1, 0)]

        def get_caps():
            mini = abs(data[0].min().min())
            maxi = abs(data[0].max().max())
            cap = round(max((mini, maxi)), 1)
            max_poss_cap = round(data[1].iloc[0].max(), 1)
            return cap, max_poss_cap
        
        def get_plot_sizes():
            nplts = [4,3]

            fig_widths = [.0001] *(nplts[1] +3)
            fig_widths[0] = driverlabels_space if driverlabels_space \
                            else config.TSH_LEFT
            if drivers_colorbar:
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
            if show_max_possible_bar:
                fig_heights[2] = config.TSH_MAX_POSS_BAR
            fig_heights[3] = config.TSH_SQUARESIZE * len(drivers)
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if targets_colorbar:
                fig_heights[4] = config.TSH_TARGETS_COLORBAR
            fig_heights[5] = config.TSH_HSPACE * (nplts[0]-1)
            fig_heights[6] = targetlabels_space if targetlabels_space \
                             else config.TSH_BOTTOM

            return nplts, fig_widths, fig_heights

        def do_plot():
            width, height = sum(fig_widths), sum(fig_heights)
            fig, axes = util._init_figure(fig_widths, fig_heights, nplts, 
                                       (config.TSH_WSPACE, config.TSH_HSPACE))
            sim, ctrl_sim = data

            if title:
                if title and type(title) is not str:
                    this_title = ('transcriptinoal similarity\nof {} and {} '
                                  .format(drivers._name, self._name))
                fig.suptitle(this_title, fontsize=config.FONTS *1.2,
                             y=1-((config.GSH_TOP/height *.25)))


            # cluster, reorder data and plot dendograms
            if cluster_targets:
                at = axes[0, 1] if show_target_dendrogram else axes[0, 0]
                order = util._heatmap_cluster(sim, 'top', at, 'columns')
                sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)
            if cluster_drivers:
                at = axes[2, 2] if show_driver_dendrogram else axes[0, 0]
                order = util._heatmap_cluster(sim, 'right', at, 'rows')
                sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order, 0)
            axes[0, 0].set_visible(False)
            
            if show_max_possible_bar:
                if reorder_to_max_possible_bar:
                    order = ctrl_sim.iloc[0].sort_values().index
                    sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)

                ctrl_lbl = drivers._ctrl if show_driverlabels else ''
                draw_cb = True if show_colorbar_legend else False
                bar_args = {'cmap': 'afmhot', 'vmin': 0, 'vmax': max_poss_cap}
                cb_lbl = 'number of markergenes' if which == 'intersect' else \
                         'base similarity\n[abs. mean eucl. distance]'
                if max_possible_bar_range is not None:
                    bar_args.update({'vmin': max_possible_bar_range[0], 
                                     'vmax': max_possible_bar_range[1]})

                util.plot_max_possible_bar(axes[1, :2], ctrl_sim,
                                            ctrl_lbl, bar_args, draw_cb, 
                                            cb_lbl, fig, width, height)

            # Y-axis setup, colorbar left
            y_ax = axes[2, 0]
            y_ax.set_ylim((-.1, sim.shape[0] +.01))
            y_ax.set_yticks(np.arange(.5, sim.shape[0]))
            ylbls = sim.index[::-1]

            if show_driverlabels:
                y_ax.tick_params(labelleft=True) 
                y_ax.set_yticklabels(ylbls, x=.5)
            if drivers_colorbar:
                prv = 0
                bar_args = config.get_plot_args('color_bar', 0)
                for col in drivers.get_colors(ylbls):
                    y_ax.bar(0, 1, color=col, bottom=prv, **bar_args)
                    prv += 1
            
            # X-axis setup, colorbar bottom
            x_ax = axes[3, 1]
            x_ax.set_xlim(0, sim.shape[1])
            x_ax.set_xticks(np.arange(.5, sim.shape[1]))
            if show_targetlabels:
                x_ax.tick_params(labelbottom=True)
                fs = targetlabel_size*config.FONTS if targetlabel_size \
                     else config.FONTS
                x_ax.set_xticklabels(sim.columns, rotation=45, ha='right',
                                     fontsize=fs, rotation_mode='anchor', y=.5)
            if targets_colorbar:
                prv = 0
                bar_args = config.get_plot_args('color_bar', 1)
                for col in self.get_colors(sim.columns):
                    x_ax.barh(0, 1, color=col, left=prv, **bar_args)
                    prv += 1
                
            # HEATMAP
            ax = axes[2, 1]
            ax.set_yticks(np.arange(0, sim.shape[0]))
            ax.set_xticks(np.arange(0, sim.shape[1]))

            hm_args = {'cmap': 'RdBu_r', 'vmin': -cap, 'vmax': cap}
            if which == 'euclid':
                if not proportional:
                    cb_lbl = 'mean change in similarity\n[mean eulc. dist.]'
                elif proportional:
                    cb_lbl = 'prop. of changed target similarity'
            elif which == 'intersect':
                if not proportional:
                    cb_lbl = ('markergene intersect\n[sum of matches(1) & '
                              'mism.(-1)]')
                elif proportional:
                    cb_lbl = 'prop. markergene intersect'
            if heatmap_range is not None:
                hm_args.update({'vmin': heatmap_range[0], 
                                'vmax': heatmap_range[1]})
            im = ax.imshow(sim.values, aspect='auto', **hm_args)
            
            # colorbar setup
            if show_colorbar_legend:
                at = (config.CB_LEFT/width, 1-config.CB_TOP/height, 
                    config.CB_WIDTH/width, config.CB_HEIGHT/height)
                cax = fig.add_axes(at)
                cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
                
                bar_ticks = [hm_args['vmin'], hm_args['vmax']]
                if which == 'intersect':
                    bar_ticks = [int(bar_ticks[0]), int(bar_ticks[1])]
                cb.set_ticks(bar_ticks)
                cb.ax.set_xticklabels(bar_ticks)
                cb.ax.set_xlabel(cb_lbl)
                cb.ax.get_xaxis().set_label_position('top')

            plt.show()
            return fig, axes, (sim, ctrl_sim)

        spacer.info('\n\n' + config.log_plot)
        check_args()
        logger.info('Processing data...')
        data = get_data()
        cap, max_poss_cap = get_caps()

        nplts, fig_widths, fig_heights = get_plot_sizes()
        logger.info('Drawing plot: {} & {}'.format(self._name, drivers._name))
        pp = util.open_file(filename)
        fig, axes, data = do_plot()
        util.save_file(fig, filename, pp)
        if pp:
            pp.close()

        return fig, axes, data
        






    def gene_similarity_heatmap(self, 
                                drivers,  
                                which,
                                differential = True,
                                proportional = False, 
                                norm_prop = True,
                                gene_number = None,
                                specific_markergenes = None,
                                custom_target_genelist = None,
                                
                                heatmap_range = None,
                                heatmap_width = None,
                                heatmap_height = None,
            
                                show_max_possible_bar = True,
                                reorder_to_max_possible_bar = False,
                                max_possible_bar_range = False,
            
                                show_sum_plot = True,
                                sum_plot_xlim = None,
                                sum_plot_central_metric = 'mean',
                                
                                cluster_genes = True,
                                show_gene_dendrogram = True,
                                genes_colorbar = None,
                                cluster_drivers = True,
                                show_driver_dendrogram = False,
                                drivers_colorbar = True,
                                
                                show_genelabels = True,
                                genelabels_space = None,
                                genelabels_size = .8,
                                bold_emphz_genes = None,
                                show_driverlabels = True,
                                driverlabels_space = None,
                                title = True, 
                                show_colorbar_legend = True,
                                filename = 'gene_similarity_hm.pdf'):

        def check_args():            
            if reorder_to_max_possible_bar:
                cluster_genes = False
            if custom_target_genelist and specific_markergenes:
                specific_markergenes = None
            assert which in ('euclid', 'intersect'), '`which` invalid '
            assert isinstance(drivers, Drivers), '`drivers` invalid'
            
            if custom_target_genelist is None:
                genes = pd.DataFrame({'name': util.annotate(self._diff_mgs.index), 
                                      'ensg': self._diff_mgs.index })
                if specific_markergenes is not None:
                    inv = list(filter(g not in genes.name, specific_markergenes))
                    specific_markergenes = list(filter(g in genes.name, 
                                                       specific_markergenes))
                    if inv:
                        logger.warning('{} are not markergenes in any of the '
                                       'targets in `{}`. These genes will not '
                                       'be included.'.format(inv, self._name))

            elif custom_target_genelist is not None:
                detec_genes = util.annotate(self._detec_genes)
                inv = list(filter(g not in detec_genes, custom_target_genelist))
                val_genes = list(filter(g in detec_genes, custom_target_genelist))
                if inv:
                    logger.warning('{} are not detected genes in any of the '
                                   'targets in `{}`. These genes will not be '
                                   'included.'.format(inv, self._name))
                    genes = pd.DataFrame({'name': val_genes, 
                                          'ensg': util.get_ensgs(val_genes)})
            return genes
                                   

        if not custom_target_genelist:
            trg = self
        else:

            cust_gl_expr = self._expr.copy().loc[genes.ensg]
            cust_gl_diff = self._diff.copy().loc[genes.ensg]
            cust_gl_diff.loc[:, :] = True
            trg = Targets(cust_gl_diff, cust_gl_expr, log=False,
                               name='custom genelist')

        sim, _, ctrl_sim, _ = trg._get_from_overlap(drivers, which, 
                                genes_diff=differential, genes_prop=proportional, 
                                norm_prop_genes=norm_prop)
        ctrl_sim = ctrl_sim.abs()

        data = {}
        def org_data(trg_sim, genes):
            trg = trg_sim.columns[0][1]
            mgt = trg_sim.columns[0][0]
            if gene_number:
                trg_sim.dropna(inplace=True)
                half = int(gene_number/2) 
                
                if which == 'euclid':
                    mins = trg_sim.min(axis=1).sort_values().index[:half]
                    maxs = trg_sim.max(axis=1).sort_values().index[-half:]
                    this_genes = genes.set_index('ensg').reindex(mins.union(maxs))
                    ts = trg_sim.reindex(this_genes.index)
                    if sum_plot_central_metric == 'mean':
                        agg_sim = ts.mean()
                    elif sum_plot_central_metric == 'median':
                        agg_sim = ts.median()
                    cs = ctrl_sim.loc[this_genes.index, (mgt, trg, drivers._ctrl)]
                
            else:
                ts = trg_sim.reindex(genes.ensg)
                if sum_plot_central_metric == 'mean':
                    agg_sim = ts.mean()
                elif sum_plot_central_metric == 'median':
                    agg_sim = ts.median()
                cs = ctrl_sim.reindex(genes.ensg)
                cs = cs.loc(1)[(mgt, trg, drivers._ctrl)]
                if not specific_markergenes_show_none:
                    ts.dropna(inplace=True)
                    cs.dropna(inplace=True)
                ord_ann_genes = genes.set_index('ensg').loc[ts.index, 'name']
            data.update({trg: {mgt: [ts.T, cs.to_frame().T, agg_sim]}})
        sim.groupby(axis=1, level=1, sort=False).apply(org_data, genes)
        data_l = [d for dat in list(data.values()) for d in list(dat.values())]

        mini = min([d.min().min() for dat in data_l for d in dat[::3]])
        maxi = max([d.max().max() for dat in data_l for d in dat[::3]])
        cap = round(max((abs(mini), abs(maxi))), 1)
        mini = min([d.min() for dat in data_l for d in dat[2::3]])
        maxi = max([d.max() for dat in data_l for d in dat[2::3]])
        agg_cap = [round(mini, 1), round(maxi, 1)]
        max_poss_cap = round(max([d.max().max() for dat in data_l 
                                 for d in dat[1::3]]), 1)
        if not differential:
            cap = max_poss_cap = max((cap, max_poss_cap))



        # =================================PLOT=================================
        pp = util.open_file(filename)
        nplts = [4, 4]
        
        fig_widths = [.0001] *(nplts[1] +3)
        fig_widths[0] = driverlabels_space if driverlabels_space else config.GSH_LEFT
        if drivers_colorbar:
            fig_widths[1] = config.GSH_DRIVERS_COLORBAR
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
        if (which == 'euclid') and show_max_possible_bar:
            fig_heights[2] = config.GSH_MAX_POSS_BAR
        fig_heights[3] = config.GSH_SQUARESIZE *len(drivers) 
        if heatmap_height:
            fig_heights[3] *= heatmap_height
        if genes_colorbar:
            fig_heights[4] = config.GSH_GENES_COLORBAR
        fig_heights[5] = config.GSH_HSPACE * (nplts[0]-1)
        fig_heights[6] = genelabels_space if genelabels_space else config.GSH_BOTTOM

        # duplicate height sizes and insert a spacer axis with the size of top
        # if self._down_mgs:
        #     nplts[0] = nplts[0] *2 +1
        #     hs = fig_heights
        #     fig_heights = hs[:-2] +[config.GSH_UP_DOWN_SPACE] +hs[1:-2] +[hs[-2:]]
        #     fig_heights[-2] = config.GSH_HSPACE * (nplts[0]-1)
        
        for t_name, dat in data.items():
            
            fig_widths[2] = config.GSH_SQUARESIZE *dat['up'][0].shape[1]
            if heatmap_width:
                fig_widths[2] *= heatmap_width
            width, height = sum(fig_widths), sum(fig_heights)
            fig, axes = util._init_figure(fig_widths, fig_heights, nplts, 
                                          (config.GSH_WSPACE, config.GSH_HSPACE))
            if trg._down_mgs:
                [ax.set_visible(False) for ax in axes[4,:]]

            if title:
                if title and type(title) is not str:
                    this_title = ('single-gene transcriptinoal similarity\nof '
                                 '{} and {} '.format(drivers._name, t_name))
                fig.suptitle(this_title, fontsize=config.FONTS *1.2,
                             y=1-((config.GSH_TOP/height *.25)))
                

            # iterate over 'up' and 'down' plot-halfs
            for (mgt, (sim, ctrl_sim, sim_agg)), r in zip(dat.items(), (0, 5)):

                # cluster and draw dendrogram
                if cluster_genes:
                    at = axes[r, 1] if show_gene_dendrogram else axes[r, 0]
                    clst = sim if not specific_markergenes_show_none else sim.dropna(1)
                    order = util._heatmap_cluster(clst, 'top', at, 'columns')
                    if specific_markergenes_show_none:
                        na_order = sim.columns[sim.isna().all(0)].tolist()
                        new_pos = list(at.get_position().bounds)
                        new_pos[2] -= len(na_order)*config.GSH_SQUARESIZE /width
                        at.set_position(new_pos)
                        order.extend(na_order)
                    sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)
                if cluster_drivers:
                    at = axes[2+r, 2] if show_driver_dendrogram else axes[r, 0]
                    clst = sim if not specific_markergenes_show_none else sim.dropna(1)
                    order = util._heatmap_cluster(clst, 'right', at, 'rows')
                    sim, sim_agg = util.align_indices([sim, sim_agg], order, 0)
                axes[r, 0].set_visible(False)

                if which == 'euclid' and show_max_possible_bar:
                    if reorder_to_max_possible_bar:
                        order = ctrl_sim.iloc[0].sort_values().index
                        sim, ctrl_sim = util.align_indices([sim, ctrl_sim], order)

                    bar_args = config.MAX_POSSIBLE_BAR_args(0)
                    if not max_possible_bar_range_from_config:
                        bar_args.update({'vmin': 0, 'vmax': max_poss_cap})
                    ctrl_lbl = drivers._ctrl if show_driverlabels else ''
                    draw_cb = True if show_colorbar_legend else False
                    cb_lbl = 'base similarity\n[abs. mean eucl. distance]'
                    util.plot_max_possible_bar(axes[1+r, :2], ctrl_sim,
                                               ctrl_lbl, bar_args, draw_cb, 
                                               cb_lbl, fig, width, height)

                # Y-axis setup, colorbar left
                y_ax = axes[2+r, 0]
                y_ax.set_ylim((-.1, sim.shape[0] +.01))
                y_ax.set_yticks(np.arange(.5, sim.shape[0]))
                ylbls = sim.index.unique(2)[::-1]
                if show_driverlabels:
                    y_ax.tick_params(labelleft=True)
                    y_ax.set_yticklabels(ylbls, x=.5)
                if drivers_colorbar:
                    prv = 0
                    bar_args = config.get_plot_args('color_bar', 0)
                    for col in drivers.get_colors(ylbls):
                        y_ax.bar(0, 1, color=col, bottom=prv, **bar_args)
                        prv += 1
                # if self._down_mgs:
                #     mgt_header = '{}-regulated genes'.format(mgt)
                #     pad = 6 if which == 'euclid' and show_max_possible_bar else 3
                #     y_ax.set_title(mgt_header, loc='right', fontweight='bold', 
                #                    pad=pad)


                # X-axis setup, annotate genes, colorbar bottom (genes)
                x_ax = axes[3+r, 1]
                x_ax.set_xlim(0, sim.shape[1])
                x_ax.set_xticks(np.arange(.5, sim.shape[1]))
                gene_ns = genes.set_index('ensg').reindex(sim.columns).name.values
                sim.columns = gene_ns
                ctrl_sim.columns = gene_ns
                if show_genelabels:
                    x_ax.tick_params(labelbottom=True)
                    fs = genelabels_size*config.FONTS if genelabels_size else config.FONTS
                    x_ax.set_xticklabels(sim.columns, rotation=45, ha='right',
                                         fontsize=fs, rotation_mode='anchor', y=.5)
                # mark_genes = self._diff_mgs.index[self._diff_mgs[(mgt, t_name)]]
                # gene_ns = genes.set_index('ensg').reindex(mark_genes).name.values
                # if specific_markergenes_show_none:
                #     [(ytl.set_fontweight('bold'), ytl.set_fontsize(config.FONTS *.9))
                #     for ytl in x_ax.get_xticklabels() if ytl.get_text() in gene_ns]

                if genes_colorbar:
                    cols = []
                    if 'default' not in genes_colorbar:
                        genes_colorbar['default'] = '#ffffff'
                    get_gene = lambda g: g if g in genes_colorbar else 'default'
                    cols = [genes_colorbar[get_gene(g)] for g in sim.columns]
                    prv = 0
                    bar_args = config.get_plot_args('color_bar', 1)
                    for col in cols:
                        x_ax.barh(0, 1, color=col, left=prv, **bar_args)
                        prv += 1

                # SUM PLOT
                if show_sum_plot:
                    ax = axes[2+r, 3]
                    ax.tick_params(labelbottom=True, bottom=True)
                    axes[3+r, 3].set_visible(False)
                    ax.set_axisbelow(True)
                    ax.xaxis.grid(alpha=0.8, linestyle='dashed')
                    
                    get_which = 0 if mgt == 'up' else 1
                    bar_args = config.get_plot_args('sum_plot', get_which)
                    xlim = bar_args.pop('xlim') 
                    if not sum_plot_xlim_from_config:
                        xlim = agg_cap
                    ax.set_xlim(xlim)
                    ax.set_ylim(0, sim_agg.shape[0])
                    yts = np.arange(.5, len(sim_agg))
                    ax.set_yticks(yts)

                    # get the correct label based on the data
                    cen = sum_plot_central_metric
                    if which == 'euclid' and differential:
                        lbl = ('{0} change in similarity\n[{0} eulc. dist.]'
                              .format(sum_plot_central_metric))
                        if proportional:
                            lbl = 'prop. of changed target similarity' 
                    elif which == 'euclid' and not differential:
                        lbl = ('abs. {} eucl. distance'
                               .format(sum_plot_central_metric))
                    elif which == 'intersect':
                        lbl = 'markergene intersect\n[sum; match=1, mismatch=-1]'
                        if proportional:
                            lbl = ('prop. markergene intersect\n'
                                   '[sum; match=1, mismatch=-1]')
                    ax.set_xlabel(lbl)

                    if differential:
                        f = lambda cen: config.colors[18] if cen < 0 else \
                                        config.color[14]
                        colors = list(map(f, sim_agg.values[::-1]))
                    else:
                        colors = config.colors[19]
                    if 'color' not in bar_args:
                        bar_args['color'] = colors
                    ax.barh(y=yts, width=sim_agg[::-1], **bar_args)

                # HEATMAP
                ax = axes[2+r, 1]
                ax.set_yticks(np.arange(0, sim.shape[0]))
                ax.set_xticks(np.arange(0, sim.shape[1]))

                # get lims, cmap and aspect ratio from config
                if which == 'euclid':
                    if differential and not proportional:
                        get_which = 0
                        cb_lbl = 'change in similarity\n[diff. eucl. dist.]'
                    if differential and proportional:
                        get_which = 1
                        cb_lbl = ('prop. of changed target similarity\n'
                                  '[prop. diff. mean eucl. dist.]')
                    if not differential:
                        get_which = 2
                        cb_lbl = ('absolute similarity\n'
                                  '[abs. mean eucl. dist.]')
                elif which == 'intersect':
                    get_which = 3
                    cb_lbl = 'similarity\n[marker gene intersects]'
                
                hm_args  = config.get_plot_args('single_gene_heatmap', get_which)
                if not heatmap_range_from_config:
                    hm_args.update({'vmax': cap, 'vmin': -cap})
                im = ax.imshow(sim.values, **hm_args)

                # colorbar setup
                if mgt == 'up' and show_colorbar_legend:
                    at = (config.CB_LEFT/width, 1-config.CB_TOP/height, 
                          config.CB_WIDTH/width, config.CB_HEIGHT/height)
                    cax = fig.add_axes(at)
                    cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
                    if which == 'intersect':
                        cb.ax.set_xticklabels(['down-regul.', 'None', 'up-regul.'])
    
                    bar_ticks = (hm_args['vmin'], hm_args['vmax'])
                    cb.set_ticks(bar_ticks)
                    cb.ax.set_xticklabels(bar_ticks)
                    cb.ax.set_xlabel(cb_lbl)
                    cb.ax.get_xaxis().set_label_position('top')

            # plt.show()
            last = False if not t_name == list(data.keys())[-1] else True
            util.save_file(fig, filename, pp, last=last)
        return fig, axes, data 
    
    


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
