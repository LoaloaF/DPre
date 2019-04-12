import pandas as pd
import sys

sys.path.append('./DPre/main')

from format_input import TARGET
from targets import Targets
from drivers import Drivers
import config





# prepare input expression table
expr = pd.read_csv('multi_inh/rsem-genes/genes_ntc_expression.tsv', sep='\t', 
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
c = Drivers(diff_genes = 'multi_inh/multi_single/deseq2/multi', 
            expression = expr, 
            ctrl = 'none', 
            name = 'epigenetic inhihbitor combinations',
            override_diff_names=True)

# a few handy features
c.set_colors(config.default_colors)
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
t.target_similarity_heatmap(c, 'euclid', proportional=False,  

              reorder_to_max_possible_bar=False, 
              cluster_targets=True, 
              cluster_drivers=False,
              
              heatmap_width=1, 
              heatmap_height=1, 
               
              heatmap_range_from_config=False,
              max_possible_bar_range_from_config=False,

              show_max_possible_bar=True, 
              show_target_dendrogram=True,
              show_driver_dendrogram=True, 

              drivers_colorbar = True, 
              targets_colorbar = True,
)

# plot the target similarity of the drivers, using the `intersect` method
t.target_similarity_heatmap(c, 'intersect', proportional=False,  

              reorder_to_max_possible_bar=False, 
              cluster_targets=True, 
              cluster_drivers=False,
              
              heatmap_width=1, 
              heatmap_height=1, 
               
              heatmap_range_from_config=False,
              max_possible_bar_range_from_config=False,

              show_max_possible_bar=True, 
              show_target_dendrogram=True,
              show_driver_dendrogram=True, 

              drivers_colorbar = True, 
              targets_colorbar = True,
)