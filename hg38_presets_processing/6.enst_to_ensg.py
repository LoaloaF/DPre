import pandas as pd
import os
import glob
os.chdir('hg38_presets_processing')


# # EXPRESSION
expr = pd.read_csv('hg38.fantom5.tsv', sep='\t', index_col='enst')
# 1. one gene is duplicated, the second occurance is droppped: ENST00000433169
expr = expr.loc(0)[~expr.index.duplicated()]
# read annotation reference 
ann = pd.read_csv('annotation/mart_export.tsv', sep='\t', index_col='Transcript stable ID')
# some transcript id's are not in the reference annotation: n = 86 from 20K genes
not_incl = expr.index.difference(ann.index)
# 2. drop those
expr = expr.reindex(expr.index.drop(not_incl))

# assign the initial ensg key
expr.index = ann.reindex(expr.index)['Gene stable ID'].values
# 3. 11 ensg keys that are not unique, drop:
expr = expr.loc(0)[~expr.index.duplicated()]

expr = expr.drop(['fantom_name', 'promoter_usage', 'name'], axis=1)
expr.to_csv('expression_base_human.tsv', sep='\t')
print('Expression index successfully changed.')

# MARKERGENES
os.chdir('som_results')
ann = pd.read_csv('../annotation/mart_export.tsv', sep='\t', index_col='Transcript stable ID')

os.makedirs('../human_marker_genes', exist_ok=True)
for file in glob.glob('*.tsv'):
    genelist = pd.read_csv(file, sep='\t', index_col=0)
    # 1.
    genelist = genelist.loc(0)[~genelist.index.duplicated()]
    not_incl = genelist.index.difference(ann.index)
    # 2.
    genelist = genelist.reindex(genelist.index.drop(not_incl))
    # assign
    genelist.index = ann.reindex(genelist.index)['Gene stable ID'].values
    # 3.
    genelist = genelist.loc(0)[~genelist.index.duplicated()]
    genelist.rename_axis('ensg').to_csv('../human_marker_genes/'+ file, sep='\t')
print('Marker genes index successfully changed.')
