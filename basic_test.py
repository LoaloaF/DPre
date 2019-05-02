import pandas as pd

from DPre import Samples, Targets, TARGET, config, color_legend

# import os
fn = 'C:/Users\LOaLoA\OneDrive\internship_SUSTech\inhibitor_project\emb\original/rsem-genes2/genes_cpm_expression.tsv'
# expr = pd.read_csv(fn, sep='\t')
# expr.set_index('ensg', inplace=True)
# # expr = expr.iloc[1::2]
# [expr.drop(c, axis=1, inplace=True) for c in expr.columns if 'mean_' not in c]
# expr = expr.reindex([
#                 'mean_embryo 2C',
#                 'mean_embryo 4C',
#                 'mean_embryo 8C',
#                'mean_embryo ICM',
#            'mean_embryo early2C',
#          'mean_embryo miioocyte',
#             'mean_embryo zygote'], 
#    axis=1)
# expr.columns = ['2C', '4C', '8C',  'ICM', 'early2C', 'miioocyte', 'zygote']
# dirr = 'C:/Users\LOaLoA\OneDrive\internship_SUSTech\inhibitor_project\emb\emb_inh\deseq2/res'
t = Targets(expression=fn, name='early embryonic fates')
# t.reorder(['miioocyte', 'zygote', 'early2C', '2C', '4C', '8C', 'ICM'])
# t.set_colors(config.colors)

# prepare input expression table

# initiate samples instance by passing the directory of deseq2 files, 
#   the pandas expression table, the name of the control, the object name used 
#   for headers and logging and finally weather the order differential names 
#   should be ignored 
# c = Samples(diff_genes = 'test_data/deseq2', 
#             expression = 'test_data/genes_ntc_expression.tsv', 
c = Samples(expression = 'test_data/genes_ntc_expression.tsv', 
            # ctrl = 'none', 
            name = 'epigenetic inhihbitor combinations',
            override_namematcher=True)

# a few handy features
# laod in a default target, don't sort the elements alphabetically, look for a 
#   preset color file to set default colors (available for `embryonic` and `all`)
# t = TARGET('mesoderm', sort=False, preset_colors=False)

# plot the target similarity of the samples, using the `euclid` method
# t.target_similarity_heatmap(c, 'euclid', pivot=True, cluster_samples=True,  cluster_targets=True)
#                           differential=True, 
#                           proportional=False,
#               reorder_to_required_effect_bar=True, 
#               display_similarity = None,
#               cluster_targets=False, 
#               cluster_samples=False,
              
#               heatmap_width=1, 
#               heatmap_height=1, 
               
#               heatmap_range=None,
#               required_effect_bar_range=None,

#               show_required_effect_bar=True, 
#               show_target_dendrogram=True,
#               show_driver_dendrogram=True, 

#               show_samples_colorbar = True, 
#               show_targets_colorbar = True,
#               filename='somth.png'mm
# )

# genes = ['Duxf3', 'Gata2', 'Zscan4b', 'Zscan4a', 'Zscan4d', 'Zscan4c', 'Zscan4e', 'Zscan4f', 'Zscan4a', 'Hsp90aa1', 'Hsp90ab1', # 2C genes
#     'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d', 'Utf1', 'Dppa2', 'Dppa3', 'Dppa4', 'Il6st', 'Foxo1', 'Prdm4', 'Zic1', 'Il6st',
#     'Epcam', 'Tdgf1', 'Cldn3', 'Cldn4', 'Ocln', 'Alpl', 'Crb3', 'Cldn11', 'Inadl','Terf1', # 'Initiation'
#     "Gdf3", "Nanog", "Zfp42", "Esrrb", "Prdm1", # 'Maturation'
#     "Utf1", "Dppa2", "Dppa3", "Dppa4", "Pecam1`", "Dnmt3l", "Klf2", "Klf5", # 'Stabilization'
#     'Nr6a1', 'Tbx21', 'Cldn11', 'Six4', "Prx", "gasdfgads", "Hist1h2bp",
#     "Tbx3", 'Nr5a1', 'Nr5a1', 'Nr5a3', 'Nr0b1', 'Gdf3', 'Nr6a1', 'Utf1', 'Lin28a', 'Sall4', 'Fbxo15', 'Snai2',
#     'Pkp1', 'Pkp2', 'Lox', 'Serpine1', 'Thbs1', 'Col6a1', 'Tgfb3', 'Tcf4', 'Ocln', 'Prdm1', 'Tgfb1', 'Klf8', 'Zeb1', 'Zeb2', 'Fn1',
#     'Pecam1', 'Nr6a1', 'Nr5a1', 'Nr5a2', 'Wnt5a', 'Wnt3a', 'Trp63', 'Trp53', 'Bmi1', 'Stat3', 'Lifr', 'Bcl3', 'Socs3', 'Bmp1ra', 'Cdc20', 'Ccnf',
#     'Stat3', 'Fgf4', 'Dnmt3b', 'Notch1', 'Grb2', 'Tcf3', 'Zic2', 'Sall1', 'Lin28a',
#     "Ncor1", "Ncor2", "Snai1", "Thy1", "Acta2", "Acta1", 'Actg1', 'Actg2', "Gapdh", "Actb", "Eef1a1",
#     "Pou5f1", "Sox2", "Nanog", "Cdh1", "Cdh2", "S100a4", "Esrrb", "Klf4", "Myc", "Lin28a",
#     "Pecam1", "Dppa5a", 'Dppa2',
#     "Rcor2", "Lin28b", "Foxo4", "Ccnd3", "Ccnb1", "Cbx1", "Dab2", "Gli1", "Tbx2", "Tgfb2",
#     ]
# cols_init = dict.fromkeys(['Duxf3', 'Gata2', 'Zscan4a', 'Zscan4b', 'Zscan4c', 'Zscan4d', 'Zscan4e', 'Zscan4f', 'Zscan4a', 'Hsp90aa1', 'Hsp90ab1', # 2C genes
#     'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d', 'Esrrb', 'Utf1', 'Dppa2', 'Dppa3', 'Dppa4', 'Il6st', 'Foxo1', 'Prdm4', 'Zic1', 'Il6st',
#     'Epcam', 'Tdgf1', 'Cldn3', 'Cldn4', 'Ocln', 'Alpl', 'Crb3', 'Cldn11', 'Inadl','Terf1'], config.colors[15])
# cols_mat = dict.fromkeys(["Gdf3", "Nanog", "Zfp42", "Esrrb", "Prdm1"], config.colors[16]) 
# cols_stab = dict.fromkeys(["Klf2", "Klf5", 'Nr6a1', 'Tbx21', 'Cldn11', 'Six4', "Prx",
#     "Tbx3", 'Nr5a1', 'Nr5a1', 'Nr5a3', 'Nr0b1', 'Gdf3', 'Nr6a1', 'Utf1', 'Lin28a', 'Sall4', 'Fbxo15', 'Snai2',
#     'Pkp1', 'Pkp2', 'Lox', 'Serpine1', 'Thbs1', 'Col6a1', 'Tgfb3', 'Tcf4', 'Ocln', 'Prdm1', 'Tgfb1', 'Klf8', 'Zeb1', 'Zeb2', 'Fn1',
#     'Pecam1', 'Nr6a1', 'Nr5a1', 'Nr5a2', 'Wnt5a', 'Wnt3a', 'Trp63', 'Trp53', 'Bmi1', 'Stat3', 'Lifr', 'Bcl3', 'Socs3', 'Bmp1ra', 'Cdc20', 'Ccnf',
#     'Stat3', 'Fgf4', 'Dnmt3b', 'Dnmt3l', 'Notch1', 'Grb2', 'Tcf3', 'Zic2', 'Sall1', 'Lin28a',
#     "Ncor1", "Ncor2", "Snai1", "Thy1", "Acta2", "Acta1", 'Actg1', 'Actg2', "Gapdh", "Actb", "Eef1a1",
#     "Pou5f1", "Sox2", "Nanog", "Cdh1", "Cdh2", "S100a4", "Esrrb", "Klf4", "Myc", "Lin28a",
#     "Pecam1", "Dppa5a", 'Dppa2', 
#     "Rcor2", "Lin28b", "Foxo4", "Ccnd3", "Ccnb1", "Cbx1", "Dab2", "Gli1", "Tbx2", "Tgfb2"], config.colors[17])
# genes_colorbar = {**cols_stab, **cols_init, **cols_mat}
t.gene_similarity_heatmap(c,  'euclid', filename='testfile.png', cluster_genes=False, show_samples_dendrogram=True, pivot=True)
#                                 differential = True,
#                                 proportional = True, 
#                                 gene_number = 20,
#                                 display_genes = 'increasing',
#                                 specific_genes = None,
#                                 custom_target_genelist = None,
                                
#                                 heatmap_range = None,
#                                 heatmap_width = 1,
#                                 heatmap_height = .5,
            
#                                 show_required_effect_bar = True,
#                                 reorder_to_required_effect_bar = True,
#                                 required_effect_bar_range = None,
            
#                                 show_sum_plot = True,
#                                 sum_plot_range = None,
#                                 sum_plot_central_metric = 'median',
#                                 show_driverlabels_sum_plot = False,
                                
#                                 cluster_genes = False,
#                                 show_gene_dendrogram = True,
#                                 genes_colorbar = True,
#                                 cluster_samples = True,
#                                 show_driver_dendrogram = True,
#                                 show_samples_colorbar = True,
                                
#                                 show_genelabels = True,
#                                 genelabel_space = None,
#                                 genelabel_size = .8,
#                                 show_driverlabels = False,
#                                 driverlabels_space_left = None,
#                                 title = True, 
#                                 show_colorbar_legend = True,
#                                 filename = 'gene_similarity_hm')


# t.ranked_similarity_barplot(c, 'euclid', 
#                             differential=False,
#                             proportional=True,
#                             rank_samples=True,
#                             display_negative=True,
                            
#                             colored_bars=True,
#                             title=True,
# )
