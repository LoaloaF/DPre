
"""

This is just to draw the SOM images for visual inspection;

"""


import numpy
from glbase3 import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm

expn = glload('trained_som.glb')

config.draw_mode = "png"
t = expn.som.tree(filename="som_tree.png", size=[6,18])

print(t)

expn.som.plot_gene_density(filename="plot_gene_density.png", topological=False)
expn.som.plot_gene_density(filename="plot_gene_density_top.png", topological=True)
expn.som.view_map(which_dim='all', pack=True, text_size=6, filename='som_test.png')
