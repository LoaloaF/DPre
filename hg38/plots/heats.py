


import matplotlib.cm as cm
from glbase3 import *
config.draw_mode = "png"

#arr = glload("../rsem-genes/genes_fc_expression.glb")
arr = glload("../hg38.fantom5.glb")
arr.log(2,.1)

z = arr.deepcopy()
z.row_Z(False)

gene_lists = {
    'H9 ESC': genelist('../som_results/H9 embryonic stem cells.tsv', format={'enst': 0}).map(genelist=arr, key='enst')['name'],
    }

gls = {}
for k in gene_lists:
    gls[k] = genelist()
    gls[k].load_list([{"name": i} for i in gene_lists[k]])

todo = {} # load prepackaged here
todo.update(gls)

for k in todo:
    mm = todo[k].map(genelist=arr, key="name")
    zz = todo[k].map(genelist=z, key="name")

    zz.heatmap(filename="heats/%s_heat.png" % k, bracket=[-2, 2],
        row_font_size=5, col_font_size=5, size=[22,22], heat_wid=0.77,
        col_cluster=True, heat_hei=0.002*len(mm),
        colbar_label="Z-score")

    mm.heatmap(filename="heats/exp_%s_heat.png" % k, bracket=[2, 10], row_norm=False,
        row_font_size=5, col_font_size=5, size=[22,22], heat_wid=0.77,
        col_cluster=True, heat_hei=0.002*len(mm),
        colbar_label="normalised expression")
