from DPre import samples
from DPre import preset_targets
from DPre import plot_color_legend
from DPre import config

# initiate samples
s = samples(expression = 'genes_cpm_expression_no_te.tsv',     # expression input data
            ctrl = 'MEF',     # name of the control in the data
            name = 'MEF to iPSC cell type conversion'    # name of the samples instance
)
s.names = ['MEF', 'OSK iCD1 Day00', 'OSK iCD1 Day01', 'OSK iCD1 Day03', 'OSK iCD1 Day05', 'OSK iCD1 Day08', 'OSK iCD1 Day10', 'ESC 2iLIF']

# initiate targets
t = preset_targets('mouse', color_legend_ncols=2)

# plot euclid target similarity
hm = t.target_similarity_heatmap(samples = s,
                                 metric = 'euclid',     # the metric to use for compu   ting similarity
                                 differential = True,
                                 heatmap_width = .067,
                                 heatmap_height = 1.45,
                                 hide_targetlabels = True,
                                 targetlabels_space = .3,
                                 hide_targets_colorbar = False,      # show the colorbar indicating target domains
                                 hide_distance_bar = False,      # hide the distance bar, for intersect n marker genes
                                 cluster_targets=True,
                                 title = False,
                                 HM_X_DENDROGRAM = .5,
                                 filename = 'target_sim_euclid.svg',     # filename for saving
                                )

# plot the absolute ranked similarity of the inital sample ans final sample
s = s.slice_elements(['MEF', 'OSK iCD1 Day10'])
t.ranked_similarity_barplot(s, 
                            metric='euclid', 
                            differential=False,
                            n_targets=6,
                            pivot=True, 
                            BP_BOTTOM=.1,
                            BP_TOP=1.5,
                            targetlabels_space=1.8,
                            title= s.names,
                            hide_base_lines = True,
                            filename='ranked.svg'
)