import pandas as pd

from DPre import Drivers, Targets, TARGET, config


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
print(c)
t.target_similarity_heatmap(c, 'euclid', proportional=False,  

              reorder_to_max_possible_bar=False, 
              cluster_targets=True, 
              cluster_drivers=False,
              
              heatmap_width=1, 
              heatmap_height=1, 
               
              heatmap_range=None,
              max_possible_bar_range=None,

              show_max_possible_bar=True, 
              show_target_dendrogram=True,
              show_driver_dendrogram=True, 

              drivers_colorbar = True, 
              targets_colorbar = True,
)

# t.gene_similarity_heatmap(c, 'euclid',

#                                 differential = True,
#                                 proportional = False, 
#                                 norm_prop = True,
#                                 gene_number = None,
#                                 specific_markergenes = None,
#                                 custom_target_genelist = None,
                                
#                                 specific_markergenes_show_none = False,
#                                 heatmap_range_from_config = False,
#                                 heatmap_width = None,
#                                 heatmap_height = None,
            
#                                 show_max_possible_bar = True,
#                                 reorder_to_max_possible_bar = False,
#                                 max_possible_bar_range_from_config = False,
            
#                                 show_sum_plot = True,
#                                 sum_plot_xlim_from_config = False,
#                                 sum_plot_central_metric = 'mean',
                                
#                                 cluster_genes = True,
#                                 show_gene_dendrogram = True,
#                                 cluster_drivers = True,
#                                 show_driver_dendrogram = False,
#                                 drivers_colorbar = True,
#                                 genes_colorbar = None,
#                                 show_colorbar_legend = True,
                                
#                                 title = True, 
#                                 show_genelabels = True,
#                                 show_driverlabels = True,
#                                 genelabels_space = None,
#                                 driverlabels_space = None,
#                                 genelabels_size = .8,
#                                 filename = 'single_gene_overlap.pdf')