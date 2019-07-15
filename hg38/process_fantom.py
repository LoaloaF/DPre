import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt

os.chdir('hg38')

expr = pd.read_csv('hg38.fantom5.tsv', sep='\t', index_col='enst')
# one gene is duplicated, the second occurance is droppped: ENST00000433169
expr.drop(expr.index[expr.index.duplicated()], axis=0, inplace=True)
ann = pd.read_csv('annotation/mart_export.tsv', sep='\t', index_col='Transcript stable ID')
# print(expr)
print(
    '\n\n'
)

# some transcript id's are not in the reference annotation: n = 86 from 20K genes
not_incl = expr.index.difference(ann.index)
# drop those
expr = expr.reindex(expr.index.drop(not_incl))

# the initial ensg key
ensg_key = ann.reindex(expr.index)['Gene stable ID'].values

# 11 ensg keys that are not unique:
dupl_ensg = ensg_key[pd.Index(ensg_key).duplicated()]
print(dupl_ensg)
# print(expr.loc(0)[dupl_ensg])
# assign ensg key
expr.index = ensg_key
# delete the 11 duplicated genes:
expr = expr.drop(dupl_ensg, axis=0)

expr = expr.drop(['fantom_name', 'promoter_usage', 'name'], axis=1)
# expr.to_csv('../preset_targets/human/expression_base.tsv', sep='\t')
# print(expr)

    
rowwise_sd = False
print(expr)
expr = np.log2(expr +1)
expr.columns = pd.MultiIndex.from_product([expr.columns, ['log2']])

m = expr.values.mean()   
s = expr.values.std() if not rowwise_sd else expr.std(1)

z_expr = expr.apply(lambda c: (c-m) /s)
z_expr.columns = pd.MultiIndex.from_product([z_expr.columns.unique(0), ['z']])
expr =  pd.concat((expr, z_expr), axis=1).reindex(expr.columns.unique(0), axis=1, level=0)

order = pd.read_pickle('../preset_targets/human/markergenes.gzip').columns.unique(1)
expr = expr.reindex(order, axis=1, level=0)
expr.rename_axis('ensg', axis=0, inplace=True)
# expr.to_csv('../preset_targets/human/expression.tsv', sep='\t')
expr.to_pickle('../preset_targets/human/expression.gzip')
print(expr)
# plt.hist(expr.xs('z',1 ,1).values.flatten(), bins=300, density=True)
# plt.show()

