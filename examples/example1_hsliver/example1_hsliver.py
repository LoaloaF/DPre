from DPre import samples
from DPre import preset_targets
from DPre import plot_color_legend
from DPre import config

# initiate samples
s = samples(expression = 'hsliver_expression.tsv',     # expression input data
            diff_genes = ('./up_genelists', 
                          './down_genelists'),    # differential gene list directories
            ctrl = 'Day00',     # name of the control in the data
            name = 'in vitro hepatic differentiated hESCs'    # name of the samples instance
)
# initiate targets
t = preset_targets('human')

# command 
# 'python ../../dpre.py -pt "human" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs"'

# reverse element order
rev_order = list(reversed(s.names))
s.reorder(rev_order)

# plot euclid target similarity
hm = t.target_similarity_heatmap(samples = s, 
                                 which = 'euclid',     # the metric to use for computing similarity
                                 heatmap_width = .08,
                                 hide_targetlabels = True,
                                 targetlabels_space=.1,
                                 hide_targets_colorbar = False,      # show the colorbar indicating target domains
                                 hide_distance_bar = True,      # hide the distance bar, for intersect n marker genes
                                 samplelabels_space = .3,
                                 title = None,
                                 pivot = True,
                                 filename = 'target_sim_euclid.png',     # filename for saving
                                )
# command
# 'python ../../dpre.py -pt "human" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs" target_sim -w "euclid" -hw 0.08 -hta -ta 0.1 -htc -hd -sa 0.3 -t "most similarity increasing genes" -pi -f "target_sim_euclid.png" '

# plot proportional intersect target similarity
hm = t.target_similarity_heatmap(samples = s, 
                                 which = 'intersect',
                                 proportional = True,     # devide score by the number of markeregenes
                                 heatmap_width = .08,
                                 hide_targetlabels = True,
                                 targetlabels_space=.1,
                                 hide_targets_colorbar = True,
                                 hide_distance_bar = False,
                                 samplelabels_space = .3,
                                 title = None,
                                 pivot = True,
                                 filename = 'target_sim_intersect.png',
                                )

# subset targets to the hepatocyte target
hep = t.slice_elements(['Hepatocyte'], name='Hepatocyte')

# create the list of genes to analyze and prepare the colorbar input
hep_bl = ['PITX2', 'KRT19', 'AFP', 'DAB2', 'FOXP4', 'SMAD3', 'FOXA1', 'HES4', 'HNF1B']
hep_cy = ['ALB', 'MXD1', 'TTR']
plot_color_legend(['heptoblast', 'heptocytes'], 
                  [config.colors[11], config.colors[12]],
                  filename='hepato_markers_legend.png'
                  )
# set up the colors for the colorbar
colorbar_cols = {**dict.fromkeys(hep_bl, config.colors[11]), 
                 **dict.fromkeys(hep_cy, config.colors[12])}

# plot the 'euclid' similarity changes of the markergenes
hm = hep.gene_similarity_heatmap(samples = s,    # the sampels to explore similarity for
                                 custom_target_genelist = hep_bl+hep_cy,    # the list of genes to be visualized
                                 heatmap_range = [-6, 6],     # range heatmap values
                                 pivot = True,
                                 heatmap_width = .7, 
                                 heatmap_height = .9, 
                                 genelabels_space = .3, 
                                 samplelabels_space_left = .4, 
                                 title = 'hepatoblast & heptocytes markergenes',
                                 show_genes_colorbar = colorbar_cols,    # the colorbar dictionary created above
                                 hide_distance_bar = True,
                                 hide_sum_plot = True,     # do not plot the summerizing plot on the right
                                 filename = 'hsliver_marker_genes.png',
                                )
# command
# 'python ../../dpre.py -pt "human" -ts "Hepatocyte" -se "./hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs" gene_sim -cu "PITX2" -cu "KRT19" -cu "AFP" -cu "DAB2" -cu "FOXP4" -cu "SMAD3" -cu "FOXA1" -cu "HES4" -cu "HNF1B" -cu "ALB" -cu "MXD1" -cu "TTR" -hr -6 -hr 6 -pi -hw 0.7 -hh 0.9 -gls 0.3 -sa 0.4 -t "hepatoblast & heptocytes markergenes" -hd -hs -f "hsliver_marker_genes.png"'

# plot the similarity changes of the most increasing genes
hm = hep.gene_similarity_heatmap(samples = s, 
    display_genes = 'increasing',      # the genes with peak positve similarity values
    gene_number = 10,        # the number of genes to extract and display
    heatmap_range = [-6, 6],
    pivot = True,
    heatmap_width = .7, 
    heatmap_height = .9, 
    genelabels_space = .3, 
    samplelabels_space_left = .4, 
    title = 'most similarity increasing genes',
    hide_distance_bar = True,
    hide_samplelabels = True,     # the label with be shared with the marker genes plot above
    hide_sum_plot = False,     # show the summary plot
    G_HM_SUMPLOT_SIZE = .65,        # adjust the size constant defined in the config module
    filename = 'hsliver_increasing_genes.png',
   )
# command 
# 'python ../../dpre.py -pt "human" -ts "Hepatocyte" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs" gene_sim -di "increasing" -ge 10 -hr -6 -hr 6 -pi -hw 0.7 -hh 0.9 -gls 0.3 -sa 0.4 -t "most similarity increasing genes" -hd -hsa -f "hsliver_increasing_genes.png"'