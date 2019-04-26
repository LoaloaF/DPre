
import os
from glbase3 import *

# extract the data;

expn = expression(filename='GSE108087_knockdowns_expr.csv', format={'name': 0}, expn='column[1:]')

# The table ab ove only has the Gene symbol, we need to add back in the ensg:
#annot = glload(os.path.expanduser('~/mm10/mm10_ensembl_v95_ensg.glb'))
#annot = annot.getColumns(['name', 'ensg'])
#annot.saveTSV('mm10_enst_annotation.tsv', key_order=['ensg', 'name']) # This is included in the repository

annot = genelist('mm10_enst_annotation.tsv', format={'force_tsv': True, 'ensg': 0, 'name':1})

print(expn)

expn = expn.sliceConditions(['shLuc.rep1', 'shLuc.rep2', 'shLuc.rep3',
    'shMcrs1.rep1', 'shMcrs1.rep2',
    'shChd4.rep1', 'shChd4.rep2'])

expn = annot.map(genelist=expn, key='name')

expn = expn.removeDuplicates('name') # R doesn't like duplicates

expn.saveTSV('raw_table.tsv')
