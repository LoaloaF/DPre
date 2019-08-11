import DPre
from DPre import preset_targets
from DPre import plot_color_legend
from DPre import config

# initiate samples 
s = DPre.samples(expression = 'hsliver_expression.tsv',     # expression input data
                 diff_genes = None,    # differential gene are loaded in the next step
                 ctrl = 'Day00',     # name of the control in the data
                 name = 'in vitro hepatic differentiated hESCs'    # name of the samples instance
)
DPre.add_diff_genes_from_z(s)

# initiate targets
t = preset_targets('human')     # initiate one of the pre-packaged targets, here the human reference

# equivalent command line (this one doesn't run because the plot argument misses)
# !python ../../dpre.py --preset_targets "human" --samples_expression "hsliver_expression.tsv" --samples_diff_genes "./up_genelists" --samples_diff_genes "./down_genelists" --control "Day00"  --samples_name "in vitro hepatic differentiated hESCs"
# !python ../../dpre.py -pt "human" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs"

# set new names
s.names = ['Day 0', 'Day 1', 'Day 2', 'Day 3', 'Day 5', 'Day 7', 'Day 9', 'Day 11', 'Day 13', 'Day 21']

# target similarity heatmap
hm = t.target_similarity_heatmap(samples = s,     # the samples to check similarity for
                                 metric = 'intersect',     # the metric to use for computing similarity
                                 filename = None,     # don't save this initial plot
                                )
hm = t.target_similarity_heatmap(samples = s,
                                 metric = 'intersect',    
                                 heatmap_width = .14, 
                                 specific_target_labels = ['Hepatocyte', 'Liver adult', 'Liver fetal'],       # only print a specifc target label instead of all
                                 targetlabels_space = .1,     # size in inches under the heatmap
#                                  filename = None,
                                )
hm = t.target_similarity_heatmap(samples = s, 
                                 metric = 'intersect',  
                                 heatmap_width = .14,
                                 specific_target_labels = ['Hepatocyte'],
                                 targetlabels_space = .4,
                                 hide_distance_bar = True,      # hide the distance bar, for intersect n marker genes
                                 cluster_targets = True,          # emphasize peak samples and targets
                                 hide_targets_dendrogram = True,  # don't show the dendrogram to save space
                                 samplelabels_space = .4,       # save some space on the left
                                 title = None,                  # remove the title to increase clarity
#                                  filename = None,             # when not passed the default filename is used
                                )

# equivalent command line
# !python ../../dpre.py -pt "human" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs" target_sim -m "intersect" -hw 0.14 -st "Hepatocyte" -ta 0.4 -hd -ct -htd -sa 0.3 -t "False"

# ## Plotting single gene transcriptional similarity
hep = t.slice_elements(['Hepatocyte'], name='Hepatocyte')
hm = hep.gene_similarity_heatmap(samples = s,              # we don't need to pass metric here, 'euclid' is the default
                                 filename = None,
                                 )

# prepare the colorbar
hep_bl = ['PITX2', 'KRT19', 'AFP', 'DAB2', 'FOXP4', 'SMAD3', 'FOXA1', 'HES4', 'HNF1B']
hep_cy = ['ALB', 'MXD1', 'TTR']
plot_color_legend(['Hepatoblast', 'Hepatocytes'], 
                  [config.colors[11], config.colors[12]],
                  filename='hepato_markers_legend.png',
                  ncolumns = 1,
                 )
# set up the colors for the colorbar
colorbar_cols = {**dict.fromkeys(hep_bl, config.colors[11]), 
                 **dict.fromkeys(hep_cy, config.colors[12])}
print(colorbar_cols)

# gene list similarity heatmap - predefined genes
hm = hep.gene_similarity_heatmap(samples = s,    # the sampels to explore similarity for
                                 custom_target_genelist = hep_bl+hep_cy,    # the list of genes to be visualized
                                 heatmap_range = [-6, 6],     # range heatmap values
                                 heatmap_width = .8, 
                                 heatmap_height = .9, 
                                 genelabels_space = .5, 
                                 samplelabels_space = .1, 
                                 title = None,
                                 hide_distance_bar = True, 
                                 show_genes_colorbar = colorbar_cols,    # the colorbar dictionary created above
                                 pivot = True,          # rotate the plot
                                 hide_samplelabels = True,         # in the final plot, shared with the previous plot
                                 hide_sum_plot = True,     # do not plot the summerizing plot on the right
                                 filename = 'hsliver_marker_genes.png',
                                )

# equivalent command line
# !python ../../dpre.py -pt "human" -ts "Hepatocyte" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs" gene_sim -cu "PITX2" -cu "KRT19" -cu "AFP" -cu "DAB2" -cu "FOXP4" -cu "SMAD3" -cu "FOXA1" -cu "HES4" -cu "HNF1B" -cu "ALB" -cu "MXD1" -cu "TTR" -hr -6 -hr 6 -pi -hw 0.8 -hh 0.9 -ge 0.3 -sa 0.1 -ge .5  -t "false" -hd -hsa -hs -f "hsliver_marker_genes.png"

# gene similarity heatmap - most increasing
hm = hep.gene_similarity_heatmap(samples = s, 
                                 display_genes = 'increasing',      # the genes with peak positve similarity values
                                 gene_number = 10,        # the number of genes to extract and display
                                 heatmap_range = [-6, 6],
                                 pivot = True,
                                 heatmap_width = .65, 
                                 heatmap_height = .9, 
                                 genelabels_space = .5, 
                                 samplelabels_space = .4, 
                                 title = 'Most similarity\nincreasing genes',
                                 hide_distance_bar = True,
                                 hide_samplelabels = True,     
                                 hide_sum_plot = True,     
                                 HM_TOP = .7,                      # the space in inches above the heatmap
                                 filename = 'hsliver_increasing_genes.png',
                                )

# equivalent command line
# !python ../../dpre.py -pt "human" -ts "Hepatocyte" -se "hsliver_expression.tsv" -sd "./up_genelists" -sd "./down_genelists" -c "Day00"  -sn "in vitro hepatic differentiated hESCs" gene_sim -gn 10 -di "increasing" -hr -6 -hr 6 -pi -hw 0.8 -hh 0.9 -ge 0.3 -sa 0.1 -ge .5  -t "Most similarity\\nincreasing genes" -hd -hs -hsa -hs -f "hsliver_marker_genes.png"