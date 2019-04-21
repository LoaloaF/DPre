

from glbase3 import *

# extract the data;

expn = expression(filename='GSE108087_knockdowns_expr.csv', format={'name': 0}, expn='column[1:]')

print(expn)

expn = expn.sliceConditions(['shLuc.rep1', 'shLuc.rep2', 'shLuc.rep3',
    'shMcrs1.rep1', 'shMcrs1.rep2',
    'shChd4.rep1', 'shChd4.rep2'])

expn = expn.removeDuplicates('name') # R doesn't like duplicates

expn.saveTSV('raw_table.tsv')
