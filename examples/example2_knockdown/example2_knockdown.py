from DPre import samples
from DPre import preset_targets
from DPre import plot_color_legend
from DPre import annotate
from DPre import config

# inititae the samples
s = samples(diff_genes = './knockdown_deseq_res', 
            expression = 'knockdown_expression.tsv',
            ctrl = 'shLuc',
            name = 'sh knockdown in mESCs',
            override_namematcher = True,    # ignore mismatches between expression and gene list names, use with care
           )
# Create the 'mouse' reference targets instance:
t = preset_targets('mouse', preset_colors=True)

# # euclid overview heatmap
t.target_similarity_heatmap(s, 
                            metric = 'euclid', 
                            hide_targetlabels = True,
                            heatmap_width = .09,
                            targetlabels_space=.8,
                            pivot=True,
                            filename =  'target_sim_euclid.png',
                            )
# command
# python ../../dpre.py -pt "mouse" -sd "./knockdown_deseq_res" -se "knockdown_expression.tsv" -c "shLuc" -sn "sh knockdown in mESCs" -so target_sim -hta -hw 0.09 -ta 0.8 -pi -f "target_sim_euclid.png"
# intersect overview heatmap
t.target_similarity_heatmap(s,
                            metric = 'intersect', 
                            hide_targetlabels = True,
                            heatmap_width = .09,
                            pivot = True,
                            targetlabels_space=.8,
                            filename = 'target_sim_intersect.png',
                            )

# draw the ranked similarity bar plot
bp = t.ranked_similarity_barplot(samples = s, 
                                 metric = 'euclid',
                                 display_negative = True,
                                 pivot = True,
                                 filename =  'ranked_sim_eucl.pdf')
# command
# python ../../dpre.py -pt "mouse" -sd "./knockdown_deseq_res" -se "knockdown_expression.tsv" -c "shLuc" -sn "sh knockdown in mESCs" -so ranked_sim -pi -din -f "ranked_sim_eucl.pdf"
bp = t.ranked_similarity_barplot(samples = s, 
                                 metric = 'intersect',
                                 display_negative = True,
                                 pivot = True,
                                 filename =  'ranked_sim_intersect.pdf')


# filter out the erythroblast targets 
eryth_names = [n for n in t.names if 'rythroblast' in n]
t_eryth = t.slice_elements(eryth_names)
# drop the Chd4 element from the samples
mcrs4 = s.slice_elements(['shMcrs1', 'shLuc'])

# run the plot without saving to get the returned plot dictionary
gene_sim_plt = t_eryth.gene_similarity_heatmap(mcrs4, 
                                               metric = 'euclid',
                                               display_genes = 'increasing', 
                                               gene_number = 80,
                                               filename = None,
                                               )

# first index the target name, 
# then element 3 of [axes, figure, data], 
# then 'up' for the marker gene type (only one)
# then element 1 of [heatmap data, distance bar data, sum plot data]
# finally the column names of this DataFrame
genelist = gene_sim_plt['Erythroblast'][2]['up'][0].columns
# annotate the ensg keys
genelist = annotate(genelist, 'mouse')

# assemble lists containing the specific hist*-groups
hist1 = []
hist2 = []
hist3 = []
hist4 = []
for gene in genelist:
    if gene.startswith('Hist1'):
        hist1.append(gene)
    elif gene.startswith('Hist2'):
        hist2.append(gene)
    elif gene.startswith('Hist3'):
        hist3.append(gene)
    elif gene.startswith('Hist4'):
        hist4.append(gene)

# create a dictionary that maps the genenames to a respective color
hist1_cols = dict.fromkeys(hist1, config.colors[10])
hist2_cols = dict.fromkeys(hist2, config.colors[11])
hist3_cols = dict.fromkeys(hist3, config.colors[12])
hist4_cols = dict.fromkeys(hist4, config.colors[13])
genes_cb = {**hist1_cols, **hist2_cols, **hist3_cols, **hist4_cols}
# plot the color legend
plot_color_legend(('Hist1', 'Hist2', 'Hist3', 'Hist4'), config.colors[10:14],
                       filename='hist_legend.png')

# plot the most increasing genes and save
data = t_eryth.gene_similarity_heatmap(mcrs4, 
                                       metric = 'euclid',
                                       display_genes = 'increasing',
                                       gene_number = 80,
                                       heatmap_width = .9,
                                       genelabels_size = .7,
                                       genelabels_space = .5,
                                       show_genes_colorbar = genes_cb,
                                       filename = 'gene_sim_incr.pdf',
                                       HM_WSPACE = .1,           # size constant defining the size between plots (left/right)
                                       plt_show = False,
                                      )
# command
#python ../../dpre.py -pt "mouse" -ts "Basophilic erythroblast" -ts "Erythroblast" -ts "Orthochromatic erythroblast" -ts "Polychromatic erythroblast" -ts "Proerythroblast" -sd "./knockdown_deseq_res" -se "knockdown_expression.tsv" -c "shLuc" -sn "sh knockdown in mESCs" -so -ss "shMcrs1" -ss "shLuc" gene_sim -di "increasing" -gn 80 -hw 0.9 -ges 0.7 -ge 0.7 -f "gene_sim_incr.pdf"
