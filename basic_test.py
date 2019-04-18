import pandas as pd

from DPre import Drivers, Targets, TARGET, config

import os

# prepare input expression table
expr = pd.read_csv(config.DPRE_PATH + 'test_data/genes_ntc_expression.tsv', sep='\t', 
                   index_col='ensg')
order = ['4C-8C DMIT', '4C BMIT', '4C DBMIT', '4C DBMI','4C DBMT',
         '8C DMVIT', '8C DMVI', '8C DMVT', '8C DVIT', 
         'DIB', 'DMB', 'DMI', 'DMVIBTCOUCGSB', 'DMV',  'DVB', 'DVI', 
         'ICM DITCO', 'ICM DTCO', 'ICM ITCO', 
         'MIB', 'MVB', 'MVI', 'VIB', 'none']
expr = expr.reindex(['mean_ESC '+ o for o in order], axis=1)
# reorder columns to match the order of the differential data (alphabetically sorted)
expr.columns = order

# initiate drivers instance by passing the directory of deseq2 files, 
#   the pandas expression table, the name of the control, the object name used 
#   for headers and logging and finally weather the order differential names 
#   should be ignored 
c = Drivers(diff_genes = config.DPRE_PATH + 'test_data/deseq2', 
            expression = expr, 
            ctrl = 'none', 
            name = 'epigenetic inhihbitor combinations',
            override_diff_names=True)

# a few handy features
c.set_colors(config.colors)
order.remove('DMVIBTCOUCGSB')
c = c.slice_elements(order)
c.reorder(['DMV',  'DVB', 'DVI', 'ICM DITCO', 'ICM DTCO', 'ICM ITCO', 
           'MIB', 'MVB', 'MVI', 'VIB', 'none', '4C-8C DMIT', '4C BMIT', 
           '4C DBMIT', '4C DBMI','4C DBMT', '8C DMVIT', '8C DMVI', 
           '8C DMVT', '8C DVIT', 'DIB', 'DMB', 'DMI'])

# laod in a default target, don't sort the elements alphabetically, look for a 
#   preset color file to set default colors (available for `embryonic` and `all`)
t = TARGET('embryonic', sort=False, colors_from_file=True)

# plot the target similarity of the drivers, using the `euclid` method
t.target_similarity_heatmap(c, 'intersect', proportional=True,  
              reorder_to_required_effect_bar=False, 
              cluster_targets=False, 
              cluster_drivers=False,
              
              heatmap_width=1, 
              heatmap_height=1, 
               
              heatmap_range=None,
              required_effect_bar_range=None,

              show_required_effect_bar=True, 
              show_target_dendrogram=True,
              show_driver_dendrogram=True, 

              show_drivers_colorbar = True, 
              show_targets_colorbar = True,
)

genes = ['Duxf3', 'Gata2', 'Zscan4b', 'Zscan4a', 'Zscan4d', 'Zscan4c', 'Zscan4e', 'Zscan4f', 'Zscan4a', 'Hsp90aa1', 'Hsp90ab1', # 2C genes
    'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d', 'Utf1', 'Dppa2', 'Dppa3', 'Dppa4', 'Il6st', 'Foxo1', 'Prdm4', 'Zic1', 'Il6st',
    'Epcam', 'Tdgf1', 'Cldn3', 'Cldn4', 'Ocln', 'Alpl', 'Crb3', 'Cldn11', 'Inadl','Terf1', # 'Initiation'
    "Gdf3", "Nanog", "Zfp42", "Esrrb", "Prdm1", # 'Maturation'
    "Utf1", "Dppa2", "Dppa3", "Dppa4", "Pecam1`", "Dnmt3l", "Klf2", "Klf5", # 'Stabilization'
    'Nr6a1', 'Tbx21', 'Cldn11', 'Six4', "Prx", "gasdfgads", "Hist1h2bp",
    "Tbx3", 'Nr5a1', 'Nr5a1', 'Nr5a3', 'Nr0b1', 'Gdf3', 'Nr6a1', 'Utf1', 'Lin28a', 'Sall4', 'Fbxo15', 'Snai2',
    'Pkp1', 'Pkp2', 'Lox', 'Serpine1', 'Thbs1', 'Col6a1', 'Tgfb3', 'Tcf4', 'Ocln', 'Prdm1', 'Tgfb1', 'Klf8', 'Zeb1', 'Zeb2', 'Fn1',
    'Pecam1', 'Nr6a1', 'Nr5a1', 'Nr5a2', 'Wnt5a', 'Wnt3a', 'Trp63', 'Trp53', 'Bmi1', 'Stat3', 'Lifr', 'Bcl3', 'Socs3', 'Bmp1ra', 'Cdc20', 'Ccnf',
    'Stat3', 'Fgf4', 'Dnmt3b', 'Notch1', 'Grb2', 'Tcf3', 'Zic2', 'Sall1', 'Lin28a',
    "Ncor1", "Ncor2", "Snai1", "Thy1", "Acta2", "Acta1", 'Actg1', 'Actg2', "Gapdh", "Actb", "Eef1a1",
    "Pou5f1", "Sox2", "Nanog", "Cdh1", "Cdh2", "S100a4", "Esrrb", "Klf4", "Myc", "Lin28a",
    "Pecam1", "Dppa5a", 'Dppa2',
    "Rcor2", "Lin28b", "Foxo4", "Ccnd3", "Ccnb1", "Cbx1", "Dab2", "Gli1", "Tbx2", "Tgfb2",
    ]
cols_init = dict.fromkeys(['Duxf3', 'Gata2', 'Zscan4a', 'Zscan4b', 'Zscan4c', 'Zscan4d', 'Zscan4e', 'Zscan4f', 'Zscan4a', 'Hsp90aa1', 'Hsp90ab1', # 2C genes
    'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d', 'Esrrb', 'Utf1', 'Dppa2', 'Dppa3', 'Dppa4', 'Il6st', 'Foxo1', 'Prdm4', 'Zic1', 'Il6st',
    'Epcam', 'Tdgf1', 'Cldn3', 'Cldn4', 'Ocln', 'Alpl', 'Crb3', 'Cldn11', 'Inadl','Terf1'], config.colors[15])
cols_mat = dict.fromkeys(["Gdf3", "Nanog", "Zfp42", "Esrrb", "Prdm1"], config.colors[16]) 
cols_stab = dict.fromkeys(["Klf2", "Klf5", 'Nr6a1', 'Tbx21', 'Cldn11', 'Six4', "Prx",
    "Tbx3", 'Nr5a1', 'Nr5a1', 'Nr5a3', 'Nr0b1', 'Gdf3', 'Nr6a1', 'Utf1', 'Lin28a', 'Sall4', 'Fbxo15', 'Snai2',
    'Pkp1', 'Pkp2', 'Lox', 'Serpine1', 'Thbs1', 'Col6a1', 'Tgfb3', 'Tcf4', 'Ocln', 'Prdm1', 'Tgfb1', 'Klf8', 'Zeb1', 'Zeb2', 'Fn1',
    'Pecam1', 'Nr6a1', 'Nr5a1', 'Nr5a2', 'Wnt5a', 'Wnt3a', 'Trp63', 'Trp53', 'Bmi1', 'Stat3', 'Lifr', 'Bcl3', 'Socs3', 'Bmp1ra', 'Cdc20', 'Ccnf',
    'Stat3', 'Fgf4', 'Dnmt3b', 'Dnmt3l', 'Notch1', 'Grb2', 'Tcf3', 'Zic2', 'Sall1', 'Lin28a',
    "Ncor1", "Ncor2", "Snai1", "Thy1", "Acta2", "Acta1", 'Actg1', 'Actg2', "Gapdh", "Actb", "Eef1a1",
    "Pou5f1", "Sox2", "Nanog", "Cdh1", "Cdh2", "S100a4", "Esrrb", "Klf4", "Myc", "Lin28a",
    "Pecam1", "Dppa5a", 'Dppa2', 
    "Rcor2", "Lin28b", "Foxo4", "Ccnd3", "Ccnb1", "Cbx1", "Dab2", "Gli1", "Tbx2", "Tgfb2"], config.colors[17])
genes_colorbar = {**cols_stab, **cols_init, **cols_mat}
t.gene_similarity_heatmap(c,  'intersect',
                                differential = True,
                                proportional = True, 
                                gene_number = 20,
                                display_genes = 'increasing',
                                specific_genes = None,
                                custom_target_genelist = None,
                                
                                heatmap_range = None,
                                heatmap_width = 1,
                                heatmap_height = .5,
            
                                show_required_effect_bar = True,
                                reorder_to_required_effect_bar = True,
                                required_effect_bar_range = None,
            
                                show_sum_plot = True,
                                sum_plot_range = None,
                                sum_plot_central_metric = 'median',
                                show_driverlabels_sum_plot = False,
                                
                                cluster_genes = False,
                                show_gene_dendrogram = True,
                                genes_colorbar = True,
                                cluster_drivers = True,
                                show_driver_dendrogram = True,
                                show_drivers_colorbar = True,
                                
                                show_genelabels = True,
                                genelabel_space = None,
                                genelabel_size = .8,
                                show_driverlabels = False,
                                driverlabels_space = None,
                                title = True, 
                                show_colorbar_legend = True,
                                filename = 'gene_similarity_hm')