
"""

Get the C/T genes and save them into som_results/

"""

from glbase3 import *

expn = glload('trained_som.glb')

# Build the lineage meaned SOMS
cond_names = expn.getConditionNames()

som_score = 0.2

# Extract all cell-type specific genes:
config.draw_mode = "png"
res = expn.som.threshold_SOM_nodes(som_score, filename="cell_lineages_thresholds.png", som_names=expn.getConditionNames(), text_size=8)

for som_name in res:
    print()
    print(som_name)
    print(res[som_name])

for som_name in res:
    res[som_name].saveTSV("som_results/%s.tsv" % som_name.replace("/", "-"))


