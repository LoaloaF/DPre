"""

SOMS

"""

import numpy
from glbase3 import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from mpl_toolkits.mplot3d import Axes3D
config.draw_mode = "png"

expn = glload('hg38.fantom5.glb')
print(expn)

labs = ['POU5F1', 'SOX2', 'NANOG', # embryonic
    'EOMES', 'HAND1', 'TFAP2A', # Other embryonic
    'SOX7', 'HNF4A', 'GATA4', 'GATA6', 'SOX17', # endoderm
    'TRP63', # Surface ectoderm
    'MYOD1', # Muscle mesoderm
    'NR2F1', 'MSX1', 'MSX2', # Neural crest
    'TAL1', 'CEBPE',# Blood mesoderm
    'SOX10','OLIG2', 'NEUROD1', 'OLIG1', 'MYT1L', 'ASCL1', # neurectoderm masters
    'DMRTB1', # Germ cell master regulators
    ]

low_expressed_threshold = 1.0

expn.draw_scatter_CV(filename='s1.scatter_CV.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-4, 2**14], vlines=[low_expressed_threshold,],
    hlines=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])

expn = expn.filter_low_expressed(2**2.0, 2) # Data does not incldue replicates here, so set this low;

expn.draw_scatter_CV(filename='s2.scatter_CV_afterLOW.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-4, 2**14], vlines=[low_expressed_threshold,],
    hlines=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])

expn = expn.filter_high_expressed(2**7, 100)

expn.draw_scatter_CV(filename='s3.scatter_CV_afterHIGH.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-4, 2**14], vlines=[low_expressed_threshold,],
    hlines=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])

expn = expn.filter_by_CV(1.5, 20)

expn.draw_scatter_CV(filename='s4.scatter_CV_afterCV.png', size=[5,3], label_fontsize=5, label_genes=labs, label_genes_key='name',
    ylims=[0,17], xlims=[2**-4, 2**14], vlines=[low_expressed_threshold,],
    hlines=[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15])

expn.abssub(5) # Without this the SOM seems to get confused by all the non-zero values
expn.log(2,.1)

expn.save('testing.glb')

expn = glload('testing.glb')
pca = expn.get_pca(feature_key_name='name')
pca.train(10)
pca.scatter(x=1, y=2, filename='scatter12.png', figsize=[18,18], label=True, adjust_labels=False)
pca.scatter(x=3, y=4, filename='scatter34.png', figsize=[18,18], label=True, adjust_labels=False)
pca.scatter(x=5, y=6, filename='scatter56.png', figsize=[18,18], label=True, adjust_labels=False)
pca.scatter(x=7, y=8, filename='scatter78.png', figsize=[18,18], label=True, adjust_labels=False)
pca.scatter(x=9, y=10, filename='scatter910.png', figsize=[18,18], label=True, adjust_labels=False)

expn.som.config(nodenames="enst", initmethod='fullpca',
    #threshold_value=sam_map.low_expressed_threshold_log2, digitize=12, # I don't use this system as I want to show the MDS on the PCA
    seed=123456, init_whiten=True,
    image_debug='debug', components=13) # >2 will invoke MDS on the PCA
expn.som.train()

# 7 works very nicely
expn.save('trained_som.glb')
