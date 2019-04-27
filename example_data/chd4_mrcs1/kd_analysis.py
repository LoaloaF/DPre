import pandas as pd
import os
import matplotlib.pyplot as plt

from DPre import Samples    # class holding the data to test
from DPre import Targets    # class holding the data to compare against
from DPre import TARGET     # function to initiate default targets
from DPre import config     # access to more settings and colors
from DPre import color_legend   # function to plot a colorlegend
from DPre import get_ensgs
from DPre import annotate

# prepare input expression table with pandas 
expr = pd.read_csv('raw_table.tsv', 
                   sep = '\t', 
                   index_col = 'ensg')
expr.drop('name', axis=1, inplace=True)
c = ['shChd4.rep1', 'shChd4.rep2']
m = ['shMcrs1.rep1','shMcrs1.rep2']
l = ['shLuc.rep1', 'shLuc.rep2', 'shLuc.rep3']
expr = pd.concat([expr[reps].mean(1) for reps in (l, m, c)], axis=1)
expr.columns =  'shChd4', 'shMcrs1', 'shLuc'

# create sample and target data 
kd_sample = Samples(diff_genes = 'deseq_res', 
                    expression = expr,
                    ctrl = 'shLuc',
                    name = 'Mcrs1 & Chd4 knockdown',
                    override_diff_names=True)
t = TARGET('all', preset_colors=True)

if not os.path.exists('dpre_output'):
    os.mkdir('dpre_output')
# generate target similarity heatmaps euclid and intersect
t.target_similarity_heatmap(kd_sample, 
                            which = 'euclid', 
                            differential = True, 
                            cluster_targets = False,
                            cluster_samples = False,
                            show_targetlabels = False,
                            heatmap_width = .16,
                            targetlabels_space= .1,
                            filename =  'dpre_output/KD_trg_sim_euclid_diff.png' 
                            )
t.target_similarity_heatmap(kd_sample, 
                            which = 'intersect', 
                            differential = True, 
                            cluster_samples = False,
                            cluster_targets = False,
                            show_targetlabels = False,
                            show_required_effect_bar = False,
                            heatmap_width = .16,
                            targetlabels_space= .1,
                            filename =  'dpre_output/KD_trg_sim_inters_diff.png' 
                            )

# generate ranked similarity plots, euclid and intersect again
t.ranked_similarity_barplot(kd_sample, 
                            which = 'euclid',
                            pivot = True,
                            filename =  'dpre_output/KD_ranked_sim_eucl_diff.pdf')
t.ranked_similarity_barplot(kd_sample, 
                            which = 'intersect',
                            pivot = True,
                            filename =  'dpre_output/KD_ranked_sim_inters_diff.pdf')

# make a colorlegend for the target colors
color_legend(list(config.default_targets_colors.values()), 
             list(config.default_targets_colors.keys()),
             filename = 'dpre_output/trg_col_legend.png')

# slice targets to erthythroblasts and samples to mcrs4
eryth_names = [n for n in t.get_names() if 'rythroblast' in n]
t_eryth = t.slice_elements(eryth_names)
mcrs4 = kd_sample.slice_elements(['shMcrs1'])

# this is a bit cheaty (copied the increased genes here after a first run)
# make a colormap for the different histone types
genes = ['Gm11168', 'Gm10721', 'Hist1h4j', 'Hist2h4', 'Hist1h3f', 'Hist1h2bn',
       'Hist2h2ac', 'Hist1h3e', 'Hist1h4f', 'Hist4h4', 'Hist1h2ac', 'Hist1h4c',
       'Hist1h2bf', 'Hist1h4k', 'Hist1h3b', 'Hist1h1e', 'Hist2h3b', 'Gm20721',
       'Hist1h2ab', 'Hist1h2be', 'Hist1h4i', 'Hist1h4a', 'Hist1h2bg',
       'Hist1h3g', 'Hist1h2bb', 'Hist1h2bl', 'Hist1h3a', 'Hist1h4h',
       'Hist1h4d', 'Hist1h2ai', 'Hist1h3i', 'Hist1h2bj', 'Hist1h2ak',
       'Hist2h3c2', 'Hist2h2aa1', 'Hist1h2bp', 'Hist1h1b', 'Hist1h1a',
       'Hist1h2br', 'Hist1h1d', 'Slc6a9', 'Hist1h3c', 'Gm10718', 'Hist2h2aa2',
       'Hist1h2af', 'Hist1h4b', 'Hist1h3h', 'Hist1h2ag', 'Shf', 'Btn1a1',
       'Ces2b', 'Hist1h2ae', 'Hist1h2bk', 'Hist1h4n', 'Cited4', 'Xkrx', 'Igf1',
       'Hist1h3d', 'Hist3h2ba', 'Dmtn', 'Dhrs11', 'Phyhip', 'Tal1', 'Pnp2',
       'Aqp1', 'Spta1', 'AA986860', 'Slc43a1', 'Slc39a8', 'Ank1', 'Otub2',
       'Ankle1', 'F630040K05Rik', 'Ccdc24', 'Hspb6', 'Khk', 'Fam221b', 'Plcd4',
       'Shc4', 'BC025920']
hist1 = []
hist2 = []
hist3 = []
hist4 = []
for gene in genes:
    if gene.startswith('Hist1'):
        hist1.append(gene)
    elif gene.startswith('Hist2'):
        hist2.append(gene)
    elif gene.startswith('Hist3'):
        hist3.append(gene)
    elif gene.startswith('Hist4'):
        hist4.append(gene)
hist1_cols = dict.fromkeys(hist1, config.colors[10])
hist2_cols = dict.fromkeys(hist2, config.colors[11])
hist3_cols = dict.fromkeys(hist3, config.colors[12])
hist4_cols = dict.fromkeys(hist4, config.colors[13])
genes_cb = {**hist1_cols, **hist2_cols, **hist3_cols, **hist4_cols}
# make a color legend
color_legend(config.colors[10:14], ('Hist1', 'Hist2', 'Hist3', 'Hist4'),
             filename='dpre_output/hist_legend.png')

# plot the most increasing genes 
data = t_eryth.gene_similarity_heatmap(mcrs4, 
                                       which = 'euclid',
                                       display_genes = 'increasing',
                                       gene_number = 80,
                                       cluster_samples = False,
                                       cluster_genes=False,
                                       heatmap_width = .9,
                                       genelabel_size = .7,
                                       genelabel_space = .3,
                                       genes_colorbar = genes_cb,
                                       filename =  'dpre_output/KD_gene_sim_incr_diff.pdf')

