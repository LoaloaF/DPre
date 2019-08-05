import pandas as pd
import numpy as np
import sys
import os
import copy
import matplotlib.pyplot as plt
from matplotlib.colors import is_color_like
from matplotlib.lines import Line2D
from scipy.spatial import distance

from DPre.main._differential import _differential
from DPre.main.samples import samples
from DPre.main._logger import spacer, logger, log_plot
import DPre.main.config as config
import DPre.main._dpre_util as util

class targets(_differential):
    """The data to compare similarity against.

    targets can hold lists of marker genes and expression data identifying a 
    collection of comparative transcriptional identities, the targets. 

    Arguments:
        marker genes (optional): Directory with deseq2 output, directories 
            with up- and (optional) down gene lists or pandas.DataFrame. 
            Defaults to None. Gene list data has an ensg key in the first column 
            or contains an 'ensg' column label. 'Up-marker genes' should list 
            genes highly expressed in targets, 'down-marker genes' are those 
            expected to be low. When passing a DataFrame, the index should 
            consist of ensg keys, the columns of a pandas.MultiIndex with 'up' 
            or 'up' & 'down' at level 0 and the element names at level 1. 
            The dtype is bool, marker genes are stored as 'True'. If None, all 
            expression values are considered as 'marker genes'.
        expression (optional): TSV expression file or pandas.DataFrame. 
            Defaults to None.
            The TSV input should have an ensg key in the first column or an ensg 
            column label. Columns `loc`, `name`, `tss_loc` and `strand` are 
            removed automatically. The data should be exclusively numerical 
            without NaN's. When passing a DataFrame, the data can be log2- 
            and z-transformed with an ensg key index and pandas.MultiIndex 
            columns with the element names at level 0, and `log2` & `z` at 
            level 1. 
        ignore_down_mgs (bool, optional): Even if found in 'marker gene' input, 
            do not use the down marker genes for the analysis. Defaults to False.
        override_namematcher (bool, optional): When both 'marker genes' and 
            'expression' passed, this overrides the element names in 
            'marker genes'. Defaults to False. When False, element names in 
            'marker genes' and 'expression' are expected to match perfectly.
        name (str, optional): Name label of the targets. Defaults to 'Targets'. 
            Used in logging and plot annotations.
        log: (bool, optional): Log the targets initiation. Defaults to True.
        
    Note:
        At least one of 'diff_genes' and 'expression' must be passed. When both 
        are passed, the inputs must have the same element order. Gene list
        data is automatically alphabetically sorted, hence the expression order
        should concur with this.
    """
    def __init__(self, markergenes=None, expression=None, name=None, 
                 ignore_down_mgs=False, override_namematcher=False, 
                 species=None, log=True):
        # call _differential __init__ method
        super().__init__(diff_genes=markergenes, expression=expression, 
                         name=name, override_namematcher=override_namematcher, 
                         log=log)
        # define if down marker genes are used, metric species and the overlap
        self._down_mgs = not ignore_down_mgs
        self._species = species 

        self._trg_sims = {}
        self._gene_sims = {}
        
        # remove down mgs from _diff if in there but not desired by the user
        if self._has_diff:
            if not self._down_mgs and 'down' in self._diff.columns.unique(0):
                self._diff.drop('down', axis=1, level=0, inplace=True)
            elif 'down' not in self._diff.columns.unique(0) and self._down_mgs:
                self._down_mgs = False
            if log:
                spacer.info('')
                n = self._diff.sum().unstack(0).reindex(self.names).to_string()
                logger.info('Number of marker genes: \n{}'.format(n))
        # inform that not passing marker genes is not recommended
        elif log:
            self._down_mgs = False
            spacer.warning('')
            logger.warning('The targets `{}` are initiated without '
                           '`marker genes`. Note that comparing against all '
                           'genes can lead to low accuracy for defining '
                           'transcriptional similarity.'.format(self.name))
        
        # _expr_mgs store
        # s a mask of _expr that holds only the marker genes
        if self._has_expr:
            expr_mgs = util._add_mg_types(self._expr.copy(), self._down_mgs)
            if self._has_diff: 
                mg_mask = lambda trg: trg.mask(~self._diff[trg.columns[0][:-1]])
                expr_mgs = expr_mgs.groupby(level=(0,1), axis=1).apply(mg_mask)
                self._expr_mgs = expr_mgs.reindex(self._mgs)
            else:
                self._expr_mgs = expr_mgs
        if log:
            spacer.info('\n\n')
        self._log_init(log)

    def __repr__(self):
        """Get a readable summary of the samples instance"""
        return ('\n=|=|= targets-instance =|=|=\nname = {};\nelements = {};\n'
                'n = {};\nmarker genes data = {};\nexpression data = {}\n'
                .format(self.name, self.names, len(self), self._has_diff, 
                        self._has_expr))
    @property
    def _mgs(self):
        """Get the genes that are at least the marker gene of one target"""
        if self._has_diff:
            return self._diff[self._diff.any(1)].index
        elif self._has_expr:
            return self._expr[self._expr.any(1)].index

    @property
    def _mg_types(self):
        """Get the marker gene types present in the targets instance"""
        return ['up', 'down'] if self._down_mgs else ['up']

    def _compute_similarity(self, samples, metric, log=True):
        """Core function computing the similarity between samples and targets
        
            This function outputs either the euclidean distance similarity for 
            metric = 'euclid' or the marker gene & diff. gene intersect for metric 
            = 'intersect'. Returns a similarity matrix that is also stored in
            the targets dictionary '_overlaps' with the id(samples) as the key.
        """
        # check marker gene detection before overlapping
        det = self.plot_detec_mgs_prop(samples, filename=None, log=log)
        det = det.reindex(self.names, level=1)
        keep = det[det.proportion >config.DROP_TARGET_DETEC_THR].index.unique(1)
        if len(keep) != len(self):
            # drop targets with too few detected genes
            dr = pd.Index(self.names).difference(keep).tolist()
            if log:
                logger.info('{} target elements dropped due to marker gene '
                            'detection proportions lower {} ''(set in '
                            'config.DROP_TARGET_DETEC_THR):\n{}'
                            .format(len(dr), config.DROP_TARGET_DETEC_THR, dr))
            self.slice_elements(keep, inplace=True, log=False)
            # self._compute_similarity(samples, metric log=log)
        if log:
            spacer.info('')
            logger.info('Computing similarity `{}` of samples `{}` and targets: '
                        '`{}` ... '.format(metric, samples.name, self.name))
        # get expression data or gene list data of samples and targets
        if metric != 'intersect':
            # get the z-data of targets marker genes and all samples genes
            trg_data = self._expr_mgs.xs('z', 1, 2)
            smp_data = samples._expr.xs('z', 1, 1)
        else:
            # get gene list data (bool) and substitute diff data with +1 for 
            # up genes, -1 for down. Samples up and down lists are merged
            smp_data = util._bool_to_int_genes(samples._diff, return_merged=True)
            diff_mgs = self._diff.reindex(self._mgs)
            trg_data = util._bool_to_int_genes(diff_mgs, trans_updown=False)
            trg_data.mask(trg_data == 0, inplace=True)
        
        def compute_trg_sim(trg_d):
            det = trg_d.index[trg_d.notna()].intersection(smp_data.index)
            trg = trg_d.reindex(det)
            smp_d = smp_data.reindex(det)
            if metric == 'cosine':
                return smp_d.apply(lambda smp: (distance.cosine(smp, trg)-1)*-1)
            elif metric == 'pearson':
                return smp_d.apply(lambda smp: (distance.correlation(smp, trg)-1)*-1)
            else:
                # additionally save per gene similarity matrix here
                if metric == 'euclid':
                    gene_sims.append(smp_d.apply(lambda smp: abs(smp-trg)))
                    return gene_sims[-1].abs().mean()
                elif metric == 'intersect':
                    # mg_inters = lambda smp_d: abs(smp_d + trg) -1
                    # def mg_inters(smp_d):
                    #     m = abs(smp_d + trg)
                    #     return m.mask(m == 2, 1).mask(m == 1, 0).mask(m == 0, -1)
                    gene_sims.append(smp_d.apply(lambda smp_d: abs(smp_d+trg) -1))
                    return gene_sims[-1].mean()
        if metric in ['euclid', 'intersect']:
            gene_sims = []
            trg_sim = trg_data.apply(compute_trg_sim)
            
            # per gene data saving
            gene_sim = pd.concat(gene_sims, axis=1, sort=False)
            idx = [*[trg_data.columns.unique(i) for i in (0,1)], smp_data.columns]
            gene_sim.columns = pd.MultiIndex.from_product(idx)
            self._gene_sims['{}-{}'.format(id(samples), metric)] = gene_sim
        else:
            trg_sim = trg_data.apply(compute_trg_sim)

        trg_sim = util._add_mgtmean(trg_sim)
        self._trg_sims['{}-{}'.format(id(samples), metric)] = trg_sim



    def _get_similarity(self, samples, metric, which_sim='target_sim', 
                          differential=False, drop_ctrl=True, 
                          inters_to_updown_not=False, log=True):
        """Access and specifically process similarity data in overlap. Returns 
            2 elements, either target similarity and ctrl target similarity or 
            per gene similarity and per gene control similarity. Option for 
            differential and absolute similarities.
        """
        # check if overlap has already been computed, if not do overlap
        try:
            key = '{}-{}'.format(id(samples), metric)
            trg_sim = self._trg_sims[key].copy()
        except KeyError:
            self._compute_similarity(samples, metric, log=log)
            trg_sim = self._trg_sims[key].copy()
        if log:
            logger.info('Selecting and processing overlap data...')

        # ensure the overlap DataFrame has the correct ordering 
        t_ord = trg_sim.columns.unique(1)
        val_t_ord = pd.Index(self.names)
        val_t_ord = val_t_ord.drop(val_t_ord.difference(t_ord))
        if t_ord.tolist() != val_t_ord.tolist():
            trg_sim = trg_sim.reindex(val_t_ord, level=1, axis=1)
        s_ord = trg_sim.index
        val_s_ord = samples.names
        if s_ord.tolist() != val_s_ord:
            trg_sim = trg_sim.reindex(val_s_ord)

        # return target similarity, i.e. the aggregated per gene similarities
        if which_sim == 'target_sim':
            sim = trg_sim
            # expression based metrics
            if metric != 'intersect':
                if samples._ctrl:
                    ctrl_sim = sim.xs(samples._ctrl, drop_level=False)
                    if drop_ctrl:
                        sim.drop(samples._ctrl, inplace=True)
                    if differential:
                        sim = sim.apply(lambda smp: smp - ctrl_sim, axis=1)
                        # return not change but the reduction in Eucl. dist.
                        if metric == 'euclid':
                            sim *= -1
                    return sim, ctrl_sim.to_frame().T
                else:
                    return sim, None
            # gene list based metrics (intersect is differential by itself)
            else:
                if samples._ctrl and drop_ctrl:
                    sim.drop(samples._ctrl, inplace=True)
                # instead of a control, the number of marker genes are returned
                n_mgs = self._diff.sum()
                n_mgs = n_mgs.append(n_mgs.groupby(level=1, axis=0).mean())
                return sim, pd.DataFrame([n_mgs.values], ['n_mgs'], sim.columns)
        
        # return gene similarity, only available for euclid and intersect metric
        elif which_sim == 'gene_sim':
            sim = self._gene_sims[key].copy()
            if metric == 'euclid':
                if samples._ctrl:
                    ctrl_sim = sim.xs(samples._ctrl, 1, 2, drop_level=False)
                    if drop_ctrl:
                        sim.drop(samples._ctrl, axis=1, level=2, inplace=True)
                    if differential:
                        c_mask = np.repeat(ctrl_sim.values, axis=1, 
                                           repeats=len(sim.columns.unique(2)))
                        sim = c_mask - sim
                    return sim, ctrl_sim
                else:
                    return sim, None
            elif metric == 'intersect':
                if samples._ctrl and drop_ctrl:
                    sim.drop(samples._ctrl, axis=1, level=2, inplace=True)
                # by default sim holds matches (1) and mismatches (-1) 
                # this option makes -1 the positive value for down-marker genes
                if inters_to_updown_not and self._down_mgs:
                    sim['down'] *= -1
                return sim, None
            
    def plot_detec_mgs_prop(self, samples, plt_show=False, 
                            filename='detec_mgs_prop.png', 
                            specific_target_labels=None, log=True):
        """Show the proportion of detected marker genes in logs and a histogram.

            Useful for adjusting the DROP_TARGET_DETEC_THR value.

            Args:
                samples (samples): The samples instance to check the detection 
                    rate for.
                plt_show (bool, optional): Directly the histogram in a new
                    window. Defaults to False.
                filename (str, optional): Filename to save the generated 
                    histogram. Defaults to None in metric case no plot is 
                    saved.
                specific_target_labels (list, optional): define a specific set 
                    of target labels to display. Defaults to None
            Returns:
                det: A DataFrame with detection values used for logging and 
                    plotting
        """
        # get proportion of detected marker genes
        if self._has_diff:
            trg_d = self._diff
        elif self._has_expr:
            cols = pd.MultiIndex.from_product((self._mg_types, self.names))
            trg_d = pd.DataFrame(True, self._expr.index, cols)
        smp_from = samples._expr if samples._has_expr else samples._diff
        smp_d = smp_from.reindex(self._mgs).notna().iloc(1)[0]
        det = trg_d.reindex(self._mgs).apply(lambda trg: trg & smp_d).sum()
        n_mgs = trg_d.sum()
        order = (det/n_mgs).sort_values().index

        # log proportion of detected marker genes
        det = pd.DataFrame({'n marker genes': n_mgs.reindex(order), 
                            'detected in samples': det.reindex(order).values, 
                            'proportion': (det/n_mgs).reindex(order).values})
        n_trgs = 10 if not len(order) <20 else int(len(order)/2)
        edges = order.droplevel(0)[:n_trgs].append(order.droplevel(0)[-n_trgs:])
        df_edges = det.loc[(slice(None), edges), :].to_string()
        if log:
            spacer.info('')
            logger.info('Detection of targets ({}) marker genes in samples data '
                        '({}): \n{}\nShown are the {} edge proportion values.'
                        .format(self.name, samples.name, df_edges, len(edges)))
        if (det['detected in samples'] == 0).all():
            trg_genes = ', '.join(self._detec_genes[:3])
            smp_genes = ', '.join(samples._expr.index[:3]) if samples._has_expr \
                        else ', '.join(samples._diff.index[:3])
            msg = ('None of the targets marker genes were detected in the '
                   'samples. This is likely due to non-matching indeces from a '
                   'species-mismatch. Targets gene index: {} ... Samples gene '
                   'index: {}. Check the input files.'
                   .format(trg_genes, smp_genes))
            logger.error(msg)
            sys.exit(1)
        # draw the plot if filename is passed, otherwise only log and return df
        if filename or plt_show:
            fig, ax = plt.subplots()
            ax.bar(np.arange(len(order)), det.proportion, edgecolor='k',
                   width=1, color=self.get_colors(order.get_level_values(1)))
            ax.hlines(config.DROP_TARGET_DETEC_THR, 0, len(self))
            ax.yaxis.grid(alpha=0.8, linestyle='dashed')
            ax.set_xlabel(self.name, fontsize=4)
            if specific_target_labels:
                xlbl = [lbl if lbl in specific_target_labels else '' 
                        for lbl in order]
                ax.set_xticks(np.arange(len(xlbl)))
                ax.set_xticklabels(xlbl, rotation=45, ha='right', 
                                rotation_mode='anchor')
            ax.set_ylabel('Proportion of detected marker genes', fontsize=4)
            tit = ('Proportion of detected {} marker genes in {}\nline = drop '
                   'threshold').format(self.name, samples.name)
            ax.set_title(tit, fontsize=6)
            if plt_show:
                plt.show()
            if filename:
                fig.savefig(fname=filename)
                logger.info('Plot saved at {}\n'
                            .format(os.path.abspath(filename)))
            plt.close()
        return det

    def target_similarity_heatmap(self, 
                                  # plot data
                                  samples, 
                                  metric = None, 
                                  differential = True,
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
                                  specific_target_labels = None,
                                  targetlabels_space = None,
                                  samplelabels_space = None,
                                  targetlabels_size = None,
                                  samplelabels_size = None,
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
                                  # others
                                  plt_show = False,
                                  filename = 'target_similarity_hm.png',
                                  **kwargs):
        """Plot the similarity of the samples with the targets in a heatmap.
        
        This gives a compact insight on transcriptional similarity with the 
        targets. Two different metrics can be picked to assess similarity: 
        'euclid' for expression inputs or 'intersect' for comparison based 
        on diff. genes/ marker genes. Differential and absolute 
        similarity values are available options for investigating the change 
        in similarity.

        Args:
            =================== Plot data options ===================
            samples (samples): the data to rate similarity for.
            metric (str, optional): the similarity metric to use. Valid 
                options are 'euclid' and 'intersect'. Defaults to None.
                'euclid' shows the mean euclidean distance towards the 
                target marker genes expression levels and requires `expression` 
                input for samples and targets. 'intersect' will show the overlap 
                between diff. sample genes and target marker genes requiring 
                gene list input. 
                When None, determine metric based on input data in targets 
                and samples. 
            differential (bool, optional): plot the differential (change in) 
                similarity between samples and targets. Defaults to True. 
                Requires a control to be passed for the 'euclid' metric. Cannot 
                be False for 'intersect'-metric.
            display_similarity (str, optional): specify the group of 
                marker genes to display similarity for. . Defaults to 'mgs mean'. 
                Valid options are 
                'mgs mean', 'mgs up', 'mgs down'. Relevent when targets are 
                initiated with 
                down-marker genes.

            =================== data ordering options ===================
            cluster_targets (bool, optional): cluster targets using the 
                euclidean distance. Defaults to False.
            cluster_samples (bool, optional): cluster samples using the 
                euclidean distance. Defaults to False.
            reorder_to_distance_bar (bool, optional): reorder the targets
                from lowest to highest base distance. Defaults to False. 
                Cannot be True when  'cluster_targets' is True aswell.
                For details, check the 'hide_distance_bar' argument. 

            =================== general visual options ===================
            pivot (bool, optional): pivot the heatmap by 90 degrees. Defaults to 
                False. Useful for fitting the heatmap on a canvas. 
            heatmap_width (float, optional): multiplier to stretch/ squeeze 
                the heatmap squares in x direction. Defaults to None. 
                Useful for very low or high number of targets. For pivot = True
                this paramter controls the height. 
            heatmap_height (float, optional): multiplier to stretch/ squeeze 
                the heatmap squares in y direction. Defaults to None. 
                Useful for very low or high number of samples. For pivot = True
                this paramter controls the width. 
            distance_bar_range (list, optional): Define the range of values
                that form the colormap for the distance bar. Defaults to
                None. The list is interpreted as, [lower_limit, upper_limit]. 
                When None, the edges are defined to cover all occuring values. 
                
            specific_target_labels (list, optional): define a specific set of
                target labels to display. Defaults to None
            targetlabels_space (float, optional): define the size in inches
                to reserve for target labels, here, the white space on the
                bottom. Defaults to None. When None, refer to the values set in 
                config.HM_BOTTOM.
            samplelabels_space (float, optional): define the size in inches
                to reserve for sample labels, here, the white space on the
                left. Defaults to None. When None, refer to the value set in 
                config.HM_LEFT.
            targetlabels_size (float, optional): multiplier for adjusting 
                target label size. Defaults to None. Useful for very high or low 
                number of targets.
            samplelabels_size (float, optional): multiplier for adjusting 
                sample label size. Defaults to None. Useful for very high or low 
                number of samples.
            title (bool, str, optional): the plot title to set. Defaults to 
                True. For True, infer the title based on plot data inputs and 
                targets/ samples name attribute. Text input will be set as 
                the general title, False hides the title.
            kwargs: modify the constants defined in config. This is used as an 
                advanced adjustment of plot element sizes and the minimum 
                required marker genes detection proportion. This heatmap may be
                adjusted by the following paramters: DROP_TARGET_DETEC_THR,
                HM_LEFT, HM_TOP, HM_RIGHT, HM_BOTTOM, HM_WSPACE, 
                HM_HSPACE, HM_Y_COLORBAR, HM_X_COLORBAR, HM_DISTANCE_BAR, 
                HM_Y_DENDROGRAM, HM_X_DENDROGRAM, HM_SQUARE_SIZE, CB_LEFT, 
                CB_LEFT_SEC, CB_TOP, CB_WIDTH, CB_HEIGHT.
            
            =================== hide/show plot elements ===================
            hide_colorbar_legend (bool, optional): Do not plot the colorbar 
                legend. Defaults to False. Applies for all colorbar_legends.
            hide_distance_bar (bool, optional): Do not plot the distance 
                bar on top of the heatmap. Defaults to False. When True, the 
                control will appear in the main heatmap. For the 'euclid' 
                metric, this bar visualizes the absolute similarity of the 
                control with the targets. For the 'intersect' metric, the number 
                of target marker genes is shown. Defaults to False.
            hide_targetlabels (bool, optional): Do not plot the target 
                labels at the bottom. Defaults to False.
            hide_targets_dendrogram (bool, optional): Do not plot the 
                targets dendrogram from clustering. Defaults to False. 
                Requires 'cluster_targets' to be True. 
            hide_targets_colorbar (bool, optional): Do not plot the targets
                colorbar on the bottom of the heatmap. Defaults to False. 
                When colors are not set for the targets using the 
                set_colors() function, colors are set to white.
            hide_samplelabels (bool, optional): Do not plot the sample 
                labels at the left. Defaults to False.
            show_samples_dendrogram (bool, optional): Plot the samples 
                dendrogram from clustering. Defaults to False. Requires 
                'cluster_samples' to be True.
            show_samples_colorbar (bool, optional): Plot the samples
                colorbar on the left of the heatmap. Defaults to False. 
                When colors are not set for the targets using the 
                set_colors() function, colors are set to white. 

            =================== others ===================
            filename (str, optional): the filename for saving the figure.
                Defaults to 'target_similarity_hm.png'. Supported filename 
                endings are .png and .pdf. If filename does not end with 
                these, the filetype is retrieved from conifg.SAVE_FORMAT.
                If None, the plot is not saved.
            plt_show (bool, optional): directly show the created plot in a new
                window. Defaults to False.
        """
        # check user input for errors and incompatibilities
        def _check_args():
            nonlocal metric
            nonlocal differential

            nonlocal cluster_targets
            nonlocal reorder_to_distance_bar
            nonlocal hide_distance_bar
            nonlocal display_similarity
            nonlocal distance_bar_range

            # check general basic input requirements
            r = util._check_args(self, samples, metric, differential, 
                                hide_distance_bar, reorder_to_distance_bar,
                                distance_bar_range, cluster_targets, 
                                display_similarity)
            metric, differential, hide_distance_bar, reorder_to_distance_bar, \
            distance_bar_range, cluster_targets, display_similarity = r
                
            config._update_consts(kwargs)
            spacer.info('')
            logger.info('Arguments passed. Getting data now ...')

        # get the specific overlap data, plot the mean of up and down mgs
        def get_data():
            sim, ctrl_sim = self._get_similarity(samples, metric, 
                                                 differential=differential,
                                                 drop_ctrl= not hide_distance_bar)
            if ctrl_sim is None:
                ctrl_sim =  pd.DataFrame(0, [0], sim.columns)
            return [sim.xs(display_similarity, 1, 0), 
                    ctrl_sim.xs(display_similarity, 1, 0)]

        # get plot lims
        def get_caps():
            # get min and max value in data, set to caps
            if heatmap_range is not None:
                low_cap, up_cap = heatmap_range
            else:
                mini = abs(data[0].min().min())
                maxi = abs(data[0].max().max())
                up_cap = round(max((mini, maxi)), 1)
                low_cap = -up_cap
            # for the distance bar, set lims to 0,max for euclid, intersect and 
            # to -1,1 for cosine, pearson
            if distance_bar_range is not None:
                low_db_cap, up_db_cap = distance_bar_range
            else: 
                if metric in ['euclid', 'intersect']:
                    up_db_cap = round(data[1].iloc[0].max(), 1)
                    low_db_cap = 0
                elif metric == 'cosine':
                    up_db_cap = 1
                    low_db_cap = -1
                elif metric == 'pearson':
                    up_db_cap = 1
                    low_db_cap = round(data[1].iloc[0].min(), 1)
            # for absolute, both distance bar and main bar sacles must be equal
            if not differential:
                if heatmap_range is not None:
                    up_db_cap = up_cap
                    low_db_cap = low_cap
                elif metric == 'euclid':
                    up_cap = up_db_cap = max((up_cap, up_db_cap))
                    low_cap = low_db_cap
                elif metric in ['cosine', 'pearson']:
                    up_cap = up_db_cap
                    low_cap = low_db_cap

            return low_cap, up_cap, low_db_cap, up_db_cap
        
        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            nplts = [4,3]
            fig_widths = [.0001] *(nplts[1] +3)
            fig_widths[0] = samplelabels_space if samplelabels_space else \
                            config.HM_LEFT

            if show_samples_colorbar:
                fig_widths[1] = config.HM_Y_COLORBAR
            fig_widths[2] = config.HM_SQUARE_SIZE * data[0].shape[1]
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
                fig_heights[2] = config.HM_DISTANCE_BAR
            fig_heights[3] = config.HM_SQUARE_SIZE * len(samples._names_noctrl)
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if not hide_targets_colorbar:
                fig_heights[4] = config.HM_X_COLORBAR
            fig_heights[5] = config.HM_HSPACE * (nplts[0]-1)
            fig_heights[6] = targetlabels_space if targetlabels_space else \
                             config.HM_BOTTOM

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
                    this_t = util._make_title(differential, metric, 
                                              samples.name, self.name)
                else:
                    this_t = title
                if not pivot:
                    fig.suptitle(this_t, y=1- (config.HM_TOP/height)*.7, 
                                 fontsize=config.FONTS)
                else:
                    axes[2, 0].set_ylabel(this_t, labelpad=10)

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
                # set order to order of sorted values in distance bar (ctrl)
                if reorder_to_distance_bar:
                    order = ctrl_sim.iloc[0].sort_values().index
                    sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
                # only draw colorbar legend if not absolute
                draw_cb = False if hide_colorbar_legend or not differential else True
                # label of the distance bar on the left
                if metric != 'intersect' and not hide_samplelabels:
                    ctrl_lbl = samples._ctrl
                else:
                    ctrl_lbl = ''
                # general metric depended labling
                bar_args = {'cmap': 'afmhot_r', 'vmin': low_db_cap, 'vmax': up_db_cap}
                if metric == 'euclid':                    
                    cb_lbl = 'Base ' + config.EUCLID_ABS
                    bar_args.update({'cmap': 'afmhot'})
                elif metric == 'cosine':
                    cb_lbl = 'Base ' + config.COSINE_ABS
                elif metric == 'pearson':
                    cb_lbl = 'Base ' + config.PEARSON_ABS
                elif metric == 'intersect':
                    cb_lbl = config.INTERSECT_DIST_BAR
                    bar_args.update({'cmap': 'afmhot'})
                util._plot_distance_bar(axes[1, :2], ctrl_sim,
                                            ctrl_lbl, bar_args, draw_cb, 
                                            cb_lbl, fig, pivot, width, height)

            # setup heatmap x,y axis, including the colorbars
            cols = self.get_colors(sim.columns) if not hide_targets_colorbar \
                   else None
            xlbl = sim.columns
            if specific_target_labels:
                xlbl = [lbl if lbl in specific_target_labels else '' for lbl in xlbl]
            util._setup_heatmap_xy('x', axes[3, 1], xlbl, pivot,
                                  hide_targetlabels, targetlabels_size, cols)
                   
            cols = samples.get_colors(sim.index[::-1]) if show_samples_colorbar \
                   else None
            util._setup_heatmap_xy('y', axes[2, 0], sim.index[::-1], pivot, 
                                   hide_samplelabels, samplelabels_size, cols)

            ax = axes[2, 1]
            ax.set_yticks(np.arange(0, sim.shape[0]))
            ax.set_xticks(np.arange(0, sim.shape[1]))
            hm_args = {'vmin': low_cap, 'vmax': up_cap}
            hm_args['cmap'] = 'RdBu_r' if differential else 'afmhot_r'
            if metric == 'euclid' and differential:
                cb_lbl = config.EUCLID_DIFF
            elif metric == 'euclid' and not differential:
                cb_lbl = config.EUCLID_ABS
                hm_args['cmap'] = 'afmhot'
            elif metric == 'cosine' and differential:
                cb_lbl = config.COSINE_DIFF
            elif metric == 'cosine' and not differential:
                cb_lbl = config.COSINE_ABS
            elif metric == 'pearson' and differential:
                cb_lbl = config.PEARSON_DIFF
            elif metric == 'pearson' and not differential:
                cb_lbl = config.PEARSON_ABS
            elif metric == 'intersect':
                cb_lbl = config.INTERSECT
            im = ax.imshow(sim.values, aspect='auto', **hm_args)
            
            # setup heatmap colorbar legend and draw
            if not hide_colorbar_legend:
                at = (config.CB_LEFT/width, 1- config.CB_TOP/height, 
                      config.CB_WIDTH/width, config.CB_HEIGHT/height)
                cax = fig.add_axes(at)
                cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
                
                bar_ticks = [hm_args['vmin'], hm_args['vmax']]
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
        low_cap, up_cap, low_db_cap, up_db_cap = get_caps()

        nplts, fig_widths, fig_heights = get_plot_sizes()
        spacer.info('')
        logger.info('Drawing...')
        if filename:
            filename, pp = util._open_file(filename)
        fig, axes, data = do_plot()
        if plt_show:
                plt.show()
        if filename:
            util._save_file(fig, filename=filename, pp=pp)
            if pp:
                pp.close()
            logger.info('Plot saved at {}/{}\n\n'
                        .format(os.path.abspath(os.curdir), filename))
        return fig, axes, data
        

    def gene_similarity_heatmap(self,
                                # plot data
                                samples,  
                                metric = None,
                                differential = True,
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
                                samplelabels_size = None,
                                samplelabels_space = None,
                                title = True, 
                                # show/ hide elements 
                                hide_colorbar_legend = False,
                                hide_distance_bar = False,
                                hide_sum_plot = False,
                                hide_genelabels = False,
                                hide_genes_dendrogram = False,
                                show_genes_colorbar = None,
                                hide_samplelabels = False,
                                show_samples_dendrogram = False,
                                show_samples_colorbar = False,
                                # others
                                filename = 'gene_similarity_hm.pdf',
                                plt_show = False,
                                **kwargs):
        """Plot the single-gene similarities of the samples with the targets 
        in an array of heatmaps.
    
        This function reveals the drivers behind target similarity shifts. 
        Besides various gene extraction options, a genelist can be passed to 
        investigate specific similarity changes. On the right of the 
        heatmap, a bar plot visualizes a summery of the gene values.
        Two different metrics can be picked to assess similarity: 
        'euclid' for expression inputs or 'intersect' for comparison based 
        on diff. genes/ marker genes. Differential and absolute gene
        similarity values are available options for investagting the change 
        in similarity. When targets were initiated with down-marker genes,
        a seperate heatmap for each marker gene type is drawn.

        Args:
            =================== Plot data options ===================
            samples (samples): the data to rate similariity for.
            metric (str, optional): the similarity metric to use. Valid 
                options are 'euclid' and 'intersect'. Defaults to None.
                'euclid' shows the mean euclidean distance towards the 
                target marker genes expression levels and requires `expression` 
                input for samples and targets. 'intersect' will show the overlap 
                between diff. sample genes and target marker genes requiring 
                gene list input. 
                When None, determine metric based on input data in targets 
                and samples. 
            differential (bool, optional): plot the differential (change in) 
                gene similarity with the targets. Defaults to True. 
                Requires a control to be passed for the 'euclid' metric. 
                Cannot be False for 'intersect'-metric. 
            display_genes (str, optional): Extract a specific set of 
                marker genes to display for each target. Defaults to 'variant'. 
                Valid options are 'variant', 'increasing', 'decreasing' when 
                differential True, and 'variant', 'distant', 'similar' for 
                differential False. To identify sample specific effects, these 
                metrics will sort according to outlaying sample values rather 
                then overall high/ low/ increasing etc. values. This is one of 
                the 3 gene selection options to choose from.
            gene_number (int, optional): The number of genes to plot for the 
                'display_genes' option. Defaults to 45. This option is 
                ignored for the two other gene selection options 
                'specific_genes' and 'custom_target_genelist'. 
            specific_genes (list, optional): Specify the marker genes to 
                display in a list of gene names. Defaults to None. A gene 
                from this list is only displayed if it is a marker gene of 
                the specifc target and detected in the samples. This option can 
                be used idependently or in combination with 'display_genes' for 
                adding specific genes of interest to the extracted ones. Genes 
                are annotated referencing enseble v.96.
            custom_target_genelist (list, optional): Specify a custom list 
                of gene names to display similarity for. Defaults to None. 
                Currently this option is only implemented for the 'euclid' 
                similarity metric. The passed genelist will be used for all 
                targets. In contrast to 'specific_genes', the genes only need to 
                be detected in the targets instead of qualifying as specific 
                target marker genes. Still, genes need to be detected in the 
                samples. Genes are annotated referencing enseble v.96.
            
            =================== data ordering options ===================
            cluster_genes (bool, optional): cluster genes using the 
                euclidean distance. Defaults to False.
            cluster_samples (bool, optional): cluster samples using the 
                euclidean distance. Defaults to False.
            reorder_to_distance_bar (bool, optional): reorder the genes
                from lowest to highest base distance. Defaults to False. 
                Cannot be True when  'cluster_genes' is True aswell.
                For details, check the 'hide_distance_bar' argument. 

            =================== general visual options ===================
            pivot (bool, optional): pivot the heatmap by 90 degrees. Defaults 
                to False. Useful for fitting the heatmap on a canvas. 
            heatmap_width (float, optional): multiplier to stretch/ squeeze 
                the heatmap squares in x direction. Defaults to None. For 
                pivot = True this paramter controls the height. 
                Useful for very low or high number of genes. 
            heatmap_height (float, optional): multiplier to stretch/ squeeze 
                the heatmap squares in y direction. Defaults to None. For 
                pivot = True this paramter controls the width.
                Useful for very low or high number of samples.
            distance_bar_range (list, optional): Define the range of values
                that form the colormap for the distance bar. Defaults to
                None. The list is interpreted as, [lower_limit, upper_limit]. 
                When None, the edges are defined to cover 90% of occuring values 
                ignoring outlayers. 
            sum_plot_range (list, optional): Define the lower- and upper 
                x-limits for the summary plot. Defaults to None. The list is 
                interpreted as, [lower_limit, upper_limit]. When None, the 
                x-limits are defined by adding 15% to the minimum and maximum
                values. 
            genelabels_space (float, optional): define the size in inches
                to reserve for gene labels, here, the white space on the
                bottom. Defaults to None. When None, refer to the values set in 
                config.HM_BOTTOM.
            samplelabels_space (float, optional): define the size in inches
                to reserve for sample labels, here, the white space on the
                left. Defaults to None. When None, refer to the value set in 
                config.HM_LEFT.
            genelabels_size (float, optional): multiplier for adjusting gene 
                label size. Defaults to None. Useful for very high or low 
                number of genes.
            samplelabels_size (float, optional): multiplier for adjusting 
                sample label size. Defaults to None. Useful for very high or low 
                number of samples.
            title (bool, str, optional): the plot title to set. Defaults to 
                True. For True, infer the title based on plot data inputs and 
                targets/ samples name attribute. Text input will be set as 
                the general title, False hides the title.
            kwargs: modify the constants defined in config. This is used as an 
                advanced adjustment of plot element sizes and the minimum 
                required marker genes detection proportion. The heatmaps may be
                adjusted by the following paramters: DROP_TARGET_DETEC_THR,
                HM_LEFT, HM_TOP, HM_RIGHT, HM_BOTTOM, HM_WSPACE, 
                HM_HSPACE, HM_Y_COLORBAR, HM_X_COLORBAR, HM_DISTANCE_BAR, 
                HM_Y_DENDROGRAM, HM_X_DENDROGRAM, HM_SQUARE_SIZE, 
                G_HM_SUMPLOT_SIZEG_HM_UPDOWN_SPACE_SIZE, CB_LEFT, 
                CB_LEFT_SEC, CB_TOP, CB_WIDTH, 
                CB_HEIGHT.
            
            =================== hide/show plot elements ===================
            hide_colorbar_legend (bool, optional): Do not plot the colorbar 
                legend. Defaults to False. Applies for all colorbar_legends.
            hide_distance_bar (bool, optional): Do not plot the distance 
                bar on top of the heatmap. Defaults to False. When True, the 
                control will appear in the main heatmap. For the 'euclid' 
                metric, this bar visualizes the absolute similarity of the 
                control with the gene. When using the intersect metric, the 
                distance bar is never displayed.
            hide_sum_plot (bool, optional): Do not generate the summary plot on
                the left visualizing the aggregated samples. Defualts to False.
                This plot relies on the sample methad as the target similarity 
                measurement, but instead of using all target marker genes, only 
                the displayed genes are used.  
            hide_genelabels (bool, optional): Do not plot the gene 
                labels at the bottom. Defaults to False.
            hide_genes_dendrogram (bool, optional): Do not plot the 
                genes dendrogram from clustering. Defaults to False. 
                Requires 'cluster_genes' to be True. 
            show_genes_colorbar (dict, bool, optional): Plot a genes colorbar on 
                the bottom of the heatmap. Defaults to None. A dictionary 
                should map gene names to colors. Mappings for genes not 
                displayed in the plot are ignored. The color for M=missing gene 
                keys is set to white. When, True and `specifc_genes` passed,
                the passed genes will be set to config.colors[1] (green). 
            hide_samplelabels (bool, optional): Do not plot the sample 
                labels at the left. Defaults to False.
            show_samples_dendrogram (bool, optional): Plot the samples 
                dendrogram from clustering. Defaults to False. Requires 
                'cluster_samples' to be True.
            show_samples_colorbar (bool, optional): Plot the samples
                colorbar on the left of the heatmap. Defaults to False. 
                When colors are not set for the targets using the 
                set_colors() function, colors are set to white. 

            =================== others ===================
            filename (str, optional): the filename for saving the figure.
                Defaults to 'gene_similarity_hm.pdf'. Supported filename 
                endings are .png and .pdf. If filename does not end with 
                these, the filetype is retrieved from conifg.SAVE_FORMAT.
                If None, the plot is not saved.
            plt_show (bool, optional): directly show each created plot in a new
                window. Defaults to False.
        """
        # check user input for errors and incompatibilities
        def _check_args():
            nonlocal metric
            nonlocal differential
            nonlocal display_genes
            nonlocal specific_genes
            nonlocal custom_target_genelist

            nonlocal hide_distance_bar
            nonlocal reorder_to_distance_bar
            nonlocal cluster_genes
            nonlocal show_genes_colorbar
            nonlocal distance_bar_range

            # check general basic input requirements
            r = util._check_args(self, samples, metric, differential, 
                                 hide_distance_bar, reorder_to_distance_bar, 
                                 distance_bar_range, cluster_genes)
            if r[0] == 'cosine':
                # default for per gene cannot be cosine, change to euclid here
                metric = 'euclid'
            _, differential, hide_distance_bar, reorder_to_distance_bar, \
            distance_bar_range, cluster_genes, _ = r

            # check main data input
            if self._species not in ['human', 'mouse']:
                logger.info('')
                logger.error('Invalid input for species: `{}`. Valid are `mouse` '
                             'and `human`. Initate targets with these species ' 
                             'to use this function.'.format(self._species))
                sys.exit(1)
            if metric == 'intersect' and not hide_distance_bar:
                hide_distance_bar = True
                logger.warning('For the intersect metric, there is no distance'
                               'bar. `hide_distance_bar` was set to True.')
            if custom_target_genelist is not None and metric == 'intersect':
                logger.error('The `custom_target_genelist` option is '
                             'currentily not implemented for the similarity '
                             'metric `intersect`. Please choose an alternative '
                             'gene selection option.')
                sys.exit(1)
            if custom_target_genelist is not None and display_genes:
                display_genes = None
                logger.info('Both `display_genes` and '
                            '`custom_target_genelist` were passed. '
                            '`display_genes` will be ignored.')
            if display_genes:
                val = ['variant', 'increasing', 'decreasing']
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
            config._update_consts(kwargs)
            
            # modify arguments for convneience 
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
                    isin = 'marker genes'
                elif custom_target_genelist is not None:
                    inp_gl = pd.Index(custom_target_genelist).drop_duplicates()
                    val_gl_ensg = self._detec_genes.intersection(samples._detec_genes)
                    isin = 'detected genes'
                    val_gl = pd.Index(util.annotate(val_gl_ensg, self._species))
                
                inv = [g for g in inp_gl if g not in val_gl]
                inp_gl = inp_gl.drop(inv) 
                if inv:
                    logger.warning('{} ({}/{}) are not {} in any of the targets'
                                   ' or are not detected in the samples. These '
                                   'genes will not be included.'.format(inv,
                                    len(inv), len(inv)+len(inp_gl), isin))
                    if len(inv) == (len(inv)+len(inp_gl)):
                        sys.exit(1)
                # update passed list
                if specific_genes is not None:
                    specific_genes = inp_gl
                elif custom_target_genelist is not None:
                    genes = util.get_ensgs(inp_gl, self._species)
                    # duplicated indicies are painful in pandas...
                    if genes.name.duplicated().any():
                        val_gl = pd.Index(genes.ensg).intersection(val_gl_ensg)
                        genes = genes.reindex(genes.index[genes.ensg.isin(val_gl)])
                        if genes.name.tolist() != inp_gl.tolist():
                            try:
                                genes = genes.set_index('name').reindex(inp_gl)
                                genes.reset_index(inplace=True)
                                genes = genes.rename({'index': 'name'}, axis=1)
                            except Exception:
                                logger.warning('Input gene order could not be'
                                               'kept because of duplicate '
                                               'gene name issues.')
            logger.info('Arguments passed. Getting data now ...')
            return genes                   

        # get the specific similarity data and pick out the genes to display
        def get_data():
            # init a new target where all genes are marker genes of all targets
            if custom_target_genelist:
                nonlocal self
                expr = self._expr.reindex(genes.ensg).copy()
                args = {'expression': expr}
                self = targets(name='custom genelist', ignore_down_mgs=True, 
                               log=False, **args)
            sim, ctrl_sim = self._get_similarity(samples, metric, 'gene_sim',
                                                 differential=differential,
                                                 drop_ctrl= not hide_distance_bar)

            # init mutable nested dict with target and markegene type keys
            data = dict((trg, dict((mgt, None) for mgt in self._mg_types))
                        for trg in self.names)
            # select genes, form the 3 data elements per-gene similarity (heatmap), 
            # ctrl_sim (distance_bar), target similarity (sumplot)
            def sel_genes(gene_sim, genes):
                mgt = gene_sim.columns[0][0]
                trg = gene_sim.columns[0][1]
                get_genes = pd.Index([])
                gene_sim.dropna(inplace=True)
                
                if display_genes:
                    # sort similarities based on passed metric, slice to gene number
                    if display_genes == 'variant':
                        idx = gene_sim.var(1).sort_values(ascending=False).index
                    elif metric == 'euclid':
                        if display_genes in ['increasing', 'distant']:
                            idx = gene_sim.max(1).sort_values(ascending=False).index
                        elif display_genes in ['decreasing', 'similar']:
                            idx = gene_sim.min(1).sort_values().index       
                    elif metric == 'intersect':
                        if display_genes.startswith('in') and mgt == 'down' or \
                        display_genes.startswith('de') and mgt == 'up':
                            asc = True
                        elif display_genes.startswith('in') and mgt == 'up' or \
                        display_genes.startswith('de') and mgt == 'down':
                            asc = False
                        idx = gene_sim.sum(1).sort_values(ascending=asc).index
                    get_genes = idx[:gene_number]
                    
                if specific_genes is not None:
                    # check if passed genelist in target marker genes add them 
                    # if not already in 
                    inp_ensg = util.get_ensgs(specific_genes, self._species).ensg
                    not_mg = filter(lambda ie: ie not in gene_sim.index, inp_ensg)
                    inv = genes.set_index('ensg').reindex(not_mg).name
                    if not inv.empty:
                        logger.info('{} not included: not marker genes of `'
                                    '{}-{}`'.format(inv.tolist(), mgt, trg))
                    add = lambda ie: not (ie in get_genes or ie in inv)
                    add_genes = pd.Index(filter(add, inp_ensg))
                    if not add_genes.empty:
                        get_genes = get_genes.append(add_genes)
                elif custom_target_genelist:
                    get_genes = genes.ensg
                
                if get_genes.empty:
                    logger.error('No genes were picked for {}-{}. Check input.'
                                 .format(mgt, trg))
                    sys.exit(1)
                # index per gene similarity to final gene list
                # per gene similarity for heatmap
                gs = gene_sim.reindex(get_genes)
                # target similarity for heatmap
                ts = gs.mean()

                # control similarity for distance bar
                if metric == 'euclid' and not hide_distance_bar:
                    cs = ctrl_sim.loc[get_genes, (mgt, trg, samples._ctrl)].to_frame().T
                else:
                    cs = None
                data[trg][mgt] = (gs.T, cs, ts)
            
            # iterate target+marker gene type
            sim.groupby(axis=1, level=(0,1), sort=False).apply(sel_genes, genes)
            return data

        # get data limits across all targets and marker gene types to plot with 
        # one consistent heatmap range 
        def get_caps():
             # unpack nested dict into the 3 plot data elements       
            data_l = [e for dat in list(data.values())
                      for d in list(dat.values()) for e in d]
            # gene sim (heatmap), ctrl sim (distance bar) target sim (sum plot)
            gs, cs, ts = [data_l[get::3] for get in (0,1,2)]
            
             # get number of genes per plot
            n_genes = [ts.shape[1] for ts in gs]
            if self._down_mgs:
                n_genes = [max(gs[i].shape[1], gs[i+1].shape[1]) 
                           for i in range(0, len(gs), 2)]

             # get sum plot limits
            if sum_plot_range is not None:
                ts_lim = sum_plot_range
            else:
                ts_min = min([sim.min() for sim in ts])
                ts_max = max([sim.max() for sim in ts])
                ts_lim = [ts_min -abs(ts_min*.15), ts_max +abs(ts_max*.15)]
                # make sure 0 is included
                if differential or True:
                    if ts_lim[0]>=0 and ts_lim[1]>=0:
                        ts_lim[ts_lim.index(min(ts_lim))] = 0
                    elif ts_lim[0]<=0 and ts_lim[1]<=0:
                        ts_lim[ts_lim.index(max(ts_lim))] = 0
            
            # get per gene heatmap range (only required for euclid)
            if metric == 'euclid':
                if heatmap_range is not None:
                    low_cap, up_cap = heatmap_range
                else:
                    mini = [sim.min().sort_values()[int(sim.shape[1]*.05)] for sim in gs]
                    maxi = [sim.max().sort_values()[int(sim.shape[1]*.95)] for sim in gs]
                    up_cap = round(max((abs(min(mini)), abs(max(maxi)))), 1)
                    low_cap = -up_cap if differential else 0
                
                # get distance bar range
                if not hide_distance_bar:
                    if distance_bar_range is not None:
                        low_db_cap, up_db_cap = distance_bar_range
                    else: 
                        up_db_cap = round(max([sim.iloc[0].sort_values()[int(sim.shape[1]*.95)]
                                            for sim in cs]), 1)
                        low_db_cap = 0

                    # make sure heatmap and distance bar ranges align
                    if not differential:
                        if heatmap_range is not None:
                            up_db_cap = up_cap
                            low_db_cap = low_cap
                        else:
                            up_cap = up_db_cap = max((up_cap, up_db_cap))
                            low_cap = low_db_cap
                else:
                    up_db_cap = low_db_cap = None
                return up_cap, low_cap, up_db_cap, low_db_cap, ts_lim, n_genes
            # for the intersect mertic, the values can only be -1, 0 and 1
            elif metric == 'intersect':
                return 1, -1, None, None, ts_lim, n_genes
        
        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            nplts = [4, 4]
            # default size of an exes is 0
            fig_widths = [.0001] *(nplts[1] +3)
            # based on parameters and config constants, set all sizes
            fig_widths[0] = samplelabels_space if samplelabels_space \
                            else config.HM_LEFT
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
            fig_widths[6] = config.HM_RIGHT

            fig_heights = [.0001] *(nplts[0] +3)
            fig_heights[0] = config.HM_TOP
            if cluster_genes and not hide_genes_dendrogram:
                fig_heights[1] = config.HM_X_DENDROGRAM
            if not hide_distance_bar:
                fig_heights[2] = config.HM_DISTANCE_BAR 
            fig_heights[3] = config.HM_SQUARE_SIZE *len(samples._names_noctrl) 
            if heatmap_height:
                fig_heights[3] *= heatmap_height
            if show_genes_colorbar:
                fig_heights[4] = config.HM_X_COLORBAR
            fig_heights[5] = config.HM_HSPACE * (nplts[0]-1)
            fig_heights[6] = genelabels_space if genelabels_space else \
                             config.HM_BOTTOM

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
                if title == True:
                    this_t = util._make_title(differential, metric,
                                              samples.name, t_name, 
                                              postf='per gene ')
                    if display_genes:
                        this_t += ' - most similarity {} genes'.format(display_genes) 
                    elif specific_genes is not None:
                        this_t += ' - list of specific marker genes'
                    elif custom_target_genelist is not None:
                        this_t += ' - custom list of genes'
                    
                elif title and isinstance(title, (list, tuple)):
                    this_t = title[i]
                else:
                    this_t = title
                if not pivot:
                    fig.suptitle(this_t, y=1- (config.HM_TOP/height)*.7,
                                 fontsize=config.FONTS)
                else:
                    row = 2 if not self._down_mgs else 7
                    axes[row, 0].set_ylabel(this_t, labelpad=10)
                
            # iterate over up and down plot-halfs
            for mgt, r in zip(self._mg_types, (0, 5)):
                sim, ctrl_sim, sim_trg = dat[mgt]

                # cluster genes/ samples and draw dendrograms
                if cluster_genes:
                    at = axes[r, 1] if not hide_genes_dendrogram else axes[r, 0]
                    order = util._heatmap_cluster(sim, 'top', at, 'columns')
                    sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
                if cluster_samples:
                    at = axes[2+r, 2] if show_samples_dendrogram else axes[r, 0]
                    order = util._heatmap_cluster(sim, 'right', at, 'rows')
                    sim, sim_trg = util._align_indices([sim, sim_trg], order, 0)
                axes[r, 0].set_visible(False)

                # draw the distance bar 
                if not hide_distance_bar and metric == 'euclid':
                    # set order to order of sorted values in distance bar (ctrl)
                    if reorder_to_distance_bar:
                        order = ctrl_sim.iloc[0].sort_values().index
                        sim, ctrl_sim = util._align_indices([sim, ctrl_sim], order)
                        
                    bar_args = {'vmin': low_db_cap, 'vmax': up_db_cap,
                                'cmap': 'afmhot'}
                    cb_lbl = config.EUCLID_ABS
                    # only draw colorbar legend if not absolute
                    if not hide_colorbar_legend and differential and mgt=='up':
                        draw_cb = True
                    else:
                        draw_cb = False
                    # label of the distance bar on the left
                    ctrl_lbl = samples._ctrl if not hide_samplelabels else ''
                    util._plot_distance_bar(axes[1+r, :2], ctrl_sim, 
                                                  ctrl_lbl, bar_args, draw_cb, 
                                                  cb_lbl, fig, pivot, width, 
                                                  height)

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
                                      hide_samplelabels, samplelabels_size, cols)
                if self._down_mgs:
                    tit = '{} marker genes'.format(mgt)
                    pad = 13 if not hide_distance_bar else 4
                    loc = 'right' if not pivot else 'left'
                    axes[2+r, 0].set_title(tit, loc=loc, fontweight='bold', 
                                           fontsize=config.FONTS, pad=pad)
                
                # draw summary plot on the right
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
                    nsmps = sim_trg.shape[0]
                    ax.set_ylim(-.1, nsmps+.1)
                    yts = np.arange(nsmps-.5, -.5, -1)
                    ax.set_yticks(yts)

                    # setup x axes
                    ax.set_xlim(ts_lim)
                    if metric == 'euclid' and differential:
                        lbl = config.EUCLID_DIFF
                    elif metric == 'euclid' and not differential:
                        lbl = config.EUCLID_ABS
                        if not hide_distance_bar:
                            base = ctrl_sim.mean(1)
                            ax.vlines(base, 0, nsmps)
                            lbl += '\n(line = base)'
                    elif metric == 'intersect':
                        lbl = config.INTERSECT
                    if not pivot:
                        if (mgt=='up' and not self._down_mgs) or \
                        (mgt=='down' and self._down_mgs):
                            ax.set_xlabel(lbl)
                    else:
                        ax.get_yaxis().set_label_position('right')
                        ax.set_ylabel(lbl, rotation=90, labelpad=5)
                        
                    # if metric == 'euclid':
                    blue = config.colors[18] 
                    red = config.colors[14]
                    cols = [red if v >0 else blue for v in sim_trg.values]
                    ax.barh(y=yts, width=sim_trg, color=cols)

                # draw heatmap
                ax = axes[2+r, 1]
                ax.set_yticks(np.arange(0, sim.shape[0]))
                ax.set_xticks(np.arange(0, sim.shape[1]))

                hm_args = {'vmin': low_cap, 'vmax': up_cap}
                if metric == 'euclid' and differential:
                        hm_args.update({'cmap': 'RdBu_r'})
                        cb_lbl = config.EUCLID_DIFF
                if metric == 'euclid' and not differential:
                        hm_args.update({'cmap': 'afmhot'})
                        cb_lbl = config.EUCLID_ABS
                elif metric == 'intersect':
                        hm_args.update({'cmap': config.RdBu_bin})
                        cb_lbl = config.INTERSECT_GENES
                im = ax.imshow(sim.values, aspect='auto', **hm_args)

                # setup heatmap colorbar legend and draw
                if mgt == 'up' and not hide_colorbar_legend:    
                    # add a new axis for the colorbar
                    at = (config.CB_LEFT/width, 1- config.CB_TOP/height, 
                          config.CB_WIDTH/width, config.CB_HEIGHT/height)
                    cax = fig.add_axes(at)
                    cb = ax.figure.colorbar(im, cax=cax, orientation='horizontal') 
                    cb.ax.set_xlabel(cb_lbl)
                    cb.ax.get_xaxis().set_label_position('top')
                    bar_ticks = [hm_args['vmin'], hm_args['vmax']]
                    cb.set_ticks(bar_ticks)                
                    if metric == 'intersect':
                        bar_ticks = ('mismatch', 'match')
                    cb.ax.set_xticklabels(bar_ticks)
                    if pivot:
                        cb.ax.tick_params(labelrotation=90)
                dat[mgt] = sim, ctrl_sim, sim_trg
            return fig, axes, dat
            
        spacer.info('\n\n' + log_plot)
        logger.info('Plot: {} & {}'.format(self.name, samples.name))
        genes = _check_args()
        data = get_data()
        up_cap, low_cap, up_db_cap, low_db_cap, ts_lim, n_genes = get_caps()

        nplts, fig_widths, fig_heights = get_plot_sizes()
        spacer.info('')
        logger.info('Drawing...')
        if filename:
            filename, pp = util._open_file(filename)
            ftype = filename[-4:]
        ret = {}
        for i, (t_name, dat) in enumerate(data.items()):
            fig, axes, dat = do_plot(i)
            spacer.info('{}/{} --- {}'.format(i+1, len(data), t_name))
            if plt_show:
                plt.show()
            ret.update({t_name: (fig, axes, dat)})
            if filename:
                this_png_fn = '{}_{}{}'.format(filename[:-4], t_name, ftype)
                util._save_file(fig, filename=this_png_fn, pp=pp)
        if filename:
            if pp:
                pp.close()
            logger.info('Plots saved at {}/{}\n\n'
                        .format(os.path.abspath(os.curdir), filename))
        return ret 

    def ranked_similarity_barplot(self,
                                  # plot data
                                  samples,
                                  metric = None, 
                                  differential = True,
                                  display_similarity = 'mgs mean',
                                  n_targets = 16,
                                  display_negative = False,
                                  # data ordering
                                  rank_samples = False,
                                  # general settings
                                  pivot = False,
                                  xlim_range = None,
                                  targetlabels_space = None,
                                  targetlabels_size = None,
                                  colored_bars = False,
                                  spines = False,
                                  title = True,
                                  # show/ hide elements
                                  hide_targetlabels = False,
                                  hide_colorbar = False,
                                  hide_base_lines = False,
                                  # others
                                  filename = 'ranked_similarity_bp.pdf',
                                  plt_show = False,
                                  **kwargs):
        """Plot the ranked similarity of the samples with the targets in a 
        barplot

            Sort the similarity values of the samples and targets to identify
            the dominating effects in the samples. Two different metrics can be 
            picked to assess similarity: 'euclid' for expression inputs or 
            'intersect' for comparison based on diff. genes/ marker genes.
            Differential and absolute similarity values are available 
            options for investagting the change in similarity.

        Args:
            =================== Plot data options ===================
            samples (samples): the data to rank similariity for.
            metric (str, optional): the similarity metric to use. Valid 
                options are 'euclid' and 'intersect'. Defaults to None.
                'euclid' shows the mean euclidean distance towards the 
                target marker genes expression levels and requires `expression` 
                input for samples and targets. 'intersect' will show the overlap 
                between diff. sample genes and target marker genes requiring 
                gene list input. 
                When None, determine metric based on input data in targets 
                and samples.             
            differential (bool, optional): plot the differential (change in) 
                similarity between samples and targets. Defaults to True. 
                Requires a control to be passed for the 'euclid' metric. Cannot 
                be False for'intersect'-metric.
            display_similarity (str, optional): specify the group of 
                marker genes to display similarity for. . Defaults to 'mgs mean'. 
                Valid options are 
                'mgs mean', 'mgs up', 'mgs down'. Relevent when targets are 
                initiated with 
                down-marker genes.
            n_targets (int, optional): the number of targets to display in each
                plot. Defaults to 16. 
            display_negative (bool, optional): display the most negative values 
                on the bottom half of the bar plot. Defaults to False. 

            =================== data ordering options ===================
            rank_samples (bool, optional): Rank the samples based on their most 
                positive value and generate the barplots in the same order. 
                Defaults to False. When False, use the default samples order.
            
            =================== general visual options ===================
            pivot (bool, optional): pivot the barplot by 90 degrees. 
                Defaults to False. Useful for fitting the barplot on a canvas. 
            xlim_range (list, optional): Define the lower- and upper 
                x-limits for the barplot. Defaults to None. The list is 
                interpreted as, [lower_limit, upper_limit]. When None, the 
                x-limits are defined by adding 15% to the minimum and maximum
                values. 
            targetlabels_space (float, optional): define the size in inches
                to reserve for target labels, here, the white space on the
                left. Defaults to None. When None, refer to the value set in 
                config.BP_LEFT.
            targetlabels_size (float, optional): multiplier for adjusting 
                target label size. Defaults to None. Useful for very high or low 
                number of targets.
            colored_bars (bool, optional): colorize negative values in blue, 
                positive ones in red. Defaults to False.
            spines (bool, optional): in addition to the bottom and left spines,
                plot the top and right ones. Defaults to False.
            title (bool, str, optional): the plot title to set. Defaults to 
                True. For True, infer the title based on plot data inputs and 
                targets/ samples name attribute. Text input will be set as 
                the general title, False hides the title.
            kwargs: modify the constants defined in config. This is used as an 
                advanced adjustment of plot element sizes and the minimum 
                required marker genes detection proportion. The barplots may be
                adjusted by the following paramters: DROP_TARGET_DETEC_THR,
                BP_LEFT, BP_TOP, BP_RIGHT, BP_BOTTOM, BP_Y_COLORBAR, 
                BP_BARSPACE, BP_BARWIDTH_SIZE.

            =================== hide/show plot elements ===================
            hide_targetlabels (bool, optional): Do not plot the target labels 
                at the left. Defaults to False.
            hide_colorbar (bool, optional): Do not plot the targets colorbar on 
                the left of the barplot. Defaults to False. 
            hide_base_lines (bool, optional): Do not show the lines marking 
                the absolute simialrity of the control, i.e. the base line. 
            
            =================== others ===================
            filename (str, optional): the filename for saving the figure.
                Defaults to 'ranked_similarity_bp.png'. Supported filename 
                endings are .png and .pdf. If filename does not end with 
                these, the filetype is retrieved from conifg.SAVE_FORMAT.
                If None, the plot is not saved.
            plt_show (bool, optional): directly show each created plot in a new
                window. Defaults to False.
            
        """
        # check user input for errors and incompatibilities around `metric` arg
        def _check_args():
            nonlocal metric
            nonlocal differential
            nonlocal n_targets
            nonlocal display_similarity

            # check general basic input requirements
            r = util._check_args(self, samples, metric, differential,  
                                 display_similarity=display_similarity)
            metric, differential, _, _, _, _, display_similarity = r
            if not n_targets or n_targets > len(self):
                n_targets = len(self)
                logger.warning('The number of targets `n_targets` was None or '
                               'greater the length of the targets. Set to all '
                               'target elements ({}).'.format(len(self)))
            config._update_consts(kwargs)            
            logger.info('Arguments passed. Getting data now ...')

        # get the aggregated overlap data for plotting, pick the targets
        def get_data():
            sim, ctrl_sim = self._get_similarity(samples, metric, 
                                                 differential=differential)
            sim = sim.xs(display_similarity, 1, 0)
            if rank_samples:
                if differential:
                    order = sim.max(1).sort_values(ascending=False).index
                else:
                    order = sim.min(1).sort_values().index
                sim = sim.reindex(order)
            
            # slice that selects the targets in the ranking
            drop = slice(int(n_targets/2), -int(n_targets/2)) if display_negative \
                   else slice(-1, n_targets-1, -1)
            
            ascend = True if metric == 'euclid' else False
            data = dict.fromkeys(sim.index, None)
            def sel_trgs(smp_row):                
                trgs = smp_row.iloc[0].sort_values(ascending=ascend)
                data[trgs.name] = trgs.drop(trgs.index[drop])
            sim.groupby(level=0).apply(sel_trgs)
            return data, ctrl_sim
        
        # get plot global limits
        def get_caps():
            if xlim_range is not None:
                return xlim_range
            else:
                maxi = max([trg_vals.max() for trg_vals in list(data.values())])
                mini = min([trg_vals.min() for trg_vals in list(data.values())])
                ext = max([abs(maxi), abs(mini)]) *.15
                lims = [mini -ext, maxi +ext]
                if lims[0]>=0 and lims[1]>=0:
                    lims[lims.index(min(lims))] = 0
                elif lims[0]<=0 and lims[1]<=0:
                    lims[lims.index(max(lims))] = 0
                return lims

        # built 2 lists with widths and heights in inches of every axes
        def get_plot_sizes():
            fig_widths = [.0001] *5
            fig_widths[0] = targetlabels_space if targetlabels_space else \
                            config.BP_LEFT   
            if not hide_colorbar:
                fig_widths[1] = config.BP_Y_COLORBAR
            fig_widths[2] = config.BP_BARSPACE
            fig_widths[3] = .04
            fig_widths[4] = config.BP_RIGHT
            
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
            if spines:
                ax.spines['right'].set_visible(True)
                ax.spines['top'].set_visible(True)

            # set plot title
            if title in ('None', 'none', 'False', 'false', 'F', 'f'):
                this_t = False
            if title:
                if title == True:
                    this_t = util._make_title(differential, metric, s_name, 
                                              self.name, pref='ranked ')
                elif title and isinstance(title, (list, tuple)):
                    this_t = title[i]
                else:
                    this_t = title
                if not pivot:
                    fig.suptitle(this_t, y=1- (config.BP_TOP/height)*.6,
                                 fontsize=config.FONTS)
                else:
                    ax.get_yaxis().set_label_position('right')
                    ax.set_ylabel(this_t, rotation=-90, labelpad=25)

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
                fs = config.FONTS*targetlabels_size if targetlabels_size else \
                     config.FONTS
                if not pivot:
                    axes[0].set_yticklabels(ylbls, fontsize=fs)
                else:
                    axes[0].set_yticklabels(ylbls, rotation=-45, ha='right', 
                                            x=-.5, rotation_mode='anchor', 
                                            fontsize=fs)
            
            # setup x axis
            xlim = lims
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

            if metric == 'euclid' and differential:
                xlbl = config.EUCLID_DIFF
            elif metric == 'euclid' and not differential:
                xlbl = config.EUCLID_ABS
            elif metric == 'cosine' and differential:
                xlbl = config.COSINE_DIFF
            elif metric == 'cosine' and not differential:
                xlbl = config.COSINE_ABS
            elif metric == 'pearson' and differential:
                xlbl = config.PEARSON_DIFF
            elif metric == 'pearson' and not differential:
                xlbl = config.PEARSON_ABS
            elif metric == 'intersect':
                xlbl = config.INTERSECT

            # for absolute euclid sim., mark the untreated base if available
            if not differential and samples._ctrl and not hide_base_lines:
                xs = ctrl_sim.loc[samples._ctrl, display_similarity]
                xs = xs.reindex(ylbls, axis=1)
                ax.vlines(xs, yts-.4, yts+.4, linewidth=.5)
                xlbl += '\n(line = base)'
            ax.set_xlabel(xlbl, labelpad=5)

            if not colored_bars:
                cols = config.colors[19]
            else:
                blue = config.colors[18] 
                red = config.colors[14]
                cols = [red if v >0 else blue for v in dat.values]
                
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
        if filename:
            filename, pp = util._open_file(filename)
            ftype = filename[-4:]
        ret = {}
        for i, (s_name, dat) in enumerate(data.items()):
            fig, axes = do_plot(i, dat)
            spacer.info('{}/{} --- {}'.format(i+1, len(data), s_name))
            if plt_show:
                plt.show()
            ret.update({s_name: (fig, axes, dat)})
            if filename:
                this_png_fn = '{}_{}.{}'.format(filename[:-4], s_name, ftype)
                util._save_file(fig, filename=this_png_fn, pp=pp)
        if filename:
            if pp:
                pp.close()
            logger.info('Plots saved at {}/{}\n\n'
                        .format(os.path.abspath(os.curdir), filename))
        return ret 