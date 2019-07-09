
"""

This script cuts the tree to divide the cells into speciifc domains;

"""


import numpy
from glbase3 import *

expn = glload('trained_som.glb') # get this one, as it has filtered many of the genes;

config.draw_mode = "png"
c = expn.tree(filename="tree.png", label_size=3, size=[3,13], cut=0.6, row_name_key='name')

# Give the clusters names: This is done kind of crude and nowhere near as precise as in the NAR.
# Here it is just to give them nice colours anyway:
clus = {
    0: 'Neurectoderm',
    1: 'Neurectoderm',
    2: 'Neurectoderm',
    3: 'Neurectoderm',
    4: 'Neurectoderm',
    5: 'Neurectoderm',
    6: 'Blood mesoderm',
    7: 'Blood mesoderm',
    8: 'Blood mesoderm',
    9: 'Blood mesoderm',
    10: 'Blood mesoderm',
    11: 'Neurectoderm',
    12: 'Mesoderm',
    13: 'Embryonic',
    14: 'Endoderm',
    15: 'Endoderm',
    16: 'Endoderm',
    17: 'Endoderm',
    18: 'Endoderm',
    19: 'Surface ectoderm',
    20: 'Endoderm',
    21: 'Mesoderm',
    22: 'Mesoderm',
    23: 'Mesoderm',
    24: 'Mesoderm'
    }


for idx, dom in enumerate(c):
    print('%s:\t%s\t%s' % (idx, len(dom), clus[idx]))

newl = []
for idx, dom in enumerate(c):
    for ct in dom['name']:
        newl.append({'name': ct, 'cluster': idx, 'domain': clus[idx]})
gl = genelist()
gl.load_list(newl)
gl.saveTSV('domains.tsv', key_order=['name', 'cluster', 'domain'])

