import pandas as pd

from DPre import Samples    # class holding the data to test
from DPre import Targets    # class holding the data to compare against
from DPre import TARGET     # function to initiate default targets
from DPre import config     # access to more settings and colors
from DPre import color_legend   # function to plot a colorlegend

# prepare input expression table with pandas 
expr = pd.read_csv('genes_cpm_expression.tsv', 
                   sep = '\t', 
                   index_col = 'ensg')
keep_columns = [
    'ivd cardio ESC',
    'ivd cardio CM',
    'ivd cardio MES',
    'ivd mneu ESCD0',
    'ivd mneu MND4',
    'ivd ra ctrl',
    'ivd ra D4'
    ]
expr = expr.reindex(keep_columns, axis=1)
# set new more readable names
new_names = [
    'ESCs (cardio)',
    'in vitro diff. cardiomyocytes',
    'mesoderm cells (partially diff.)',
    'ESCs (mneu)',
    'in vitro diff. motor neurons',
    'ESCs (ra)',
    'ESCs +retinoic acid'
    ]
expr.columns = new_names


# initatie the target to compare against, the default target `all mouse lineages`
t = TARGET('all', sort=False, preset_colors=True)

# iniatite multiple sample data instances because of different controls
# cardiomyocytes
cm_expr = expr[['ESCs (cardio)', 'in vitro diff. cardiomyocytes', 'mesoderm cells (partially diff.)',]]
cm_sample = Samples(expression = cm_expr, 
                    ctrl = 'ESCs (cardio)',
                    name = 'in vitro differentiated ESCs')

# motor neurons
mn_expr = expr[['ESCs (mneu)', 'in vitro diff. motor neurons',]]
mn_sample = Samples(expression = mn_expr, 
                    ctrl = 'ESCs (mneu)',
                    name = 'in vitro differentiated ESCs')

# retinoic acid treated
ra_expr = expr[['ESCs (ra)', 'ESCs +retinoic acid']]
ra_sample = Samples(expression = ra_expr, 
                    ctrl = 'ESCs (ra)',
                    name = 'in vitro differentiated ESCs')

# draw the ranked similarity plot iteratively
samples = (cm_sample, mn_sample, ra_sample)
for i in range(3):
    t.ranked_similarity_barplot(samples = samples[i], 
                                which = 'euclid', 
                                differential = True,
                                n_targets = 20, 
                                show_negative = True,
                                filename = 'IVD_ranked_sim_diff.png')

# make a Samples instance with the 3 samples above to produce one heatmap
all_expr = expr[['in vitro diff. cardiomyocytes', 'in vitro diff. motor neurons', 'ESCs +retinoic acid']]
all_sample = Samples(expression = all_expr, 
                     name = 'in vitro differentiated ESCs')
t.target_similarity_heatmap(samples = all_sample, 
                            which = 'euclid', 
                            differential = False, 
                            show_driver_dendrogram = False, 
                            cluster_targets = False, 
                            cluster_samples = False,
                            heatmap_width=.2, 
                            show_targetlabels = False, 
                            targetlabels_space = .1,
                            filename = 'IVD_sim_hm_nodiff.png')

# make acolorlegend for the target colors
color_legend(list(config.default_targets_colors.values()), 
             list(config.default_targets_colors.keys()),
             filename = 'trg_legend.png')
