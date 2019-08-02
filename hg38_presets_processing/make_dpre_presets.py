"""script for creating the preset targets 'mouse' and 'human' and its sub-domains.
Genelist files generated with SOMs are read in as markergenes. Then the expression 
data is loaded and log2 z-transformed. Finally, the datasets are split in its 
domains and saved.
"""

import numpy as np
import pandas as pd
import os
import glob
import matplotlib.pyplot as plt
from DPre.main._dpre_util import _add_log2_z
from DPre.main._format_input import _format_diff_genes
from DPre.main._format_input import _format_expr
from DPre.main.config import preset_targets_colors
os.chdir('hg38_presets_processing')

species = 'human', 'mouse'
os.makedirs('../preset_targets_new', exist_ok=True)
for sp in species:
    markergenes = _format_diff_genes(sp+'_marker_genes')
    expr = _format_expr('expression_base_'+sp+'.tsv', 'targets')

    domains_map = pd.read_csv('domains_{}.tsv'.format(sp), sep='\t', index_col='name')
    domains = np.unique(domains_map.domain.values)
    # order for whole mouse/ human that matches order in specific domains
    order = pd.Index([])
    colors = []
    for dom in domains:
        dom_dir = '../preset_targets_new/{} {}'.format(sp[0], dom.lower())
    
        os.makedirs(dom_dir, exist_ok=True)
        cts = domains_map[domains_map.domain == dom].index.sort_values()
        if dom == 'Embryonic' and sp == 'mouse':
            # custom color map here
            emb_colors = pd.read_csv('mouse_embryonic_order.tsv', sep='\t', index_col=0)
            # order = development timeline
            cts = emb_colors.index.rename('name')
            pd.DataFrame(None, index=cts).to_csv('{}/{}.tsv'.format(dom_dir, dom.lower()), sep='\t')
            emb_colors.to_csv('{}/colors.tsv'.format(dom_dir), sep='\t')

        # assemble color list for whole mouse/human
        order = order.append(cts)
        color = preset_targets_colors['{} {}'.format(sp[0], dom.lower())]
        colors.append(pd.Series(color, index=cts, name='color'))
        # write final files for domains
        expr.reindex(cts, axis=1, level=0).to_pickle('{}/expression.gzip'.format(dom_dir))
        markergenes.reindex(cts, axis=1, level=1).to_pickle('{}/markergenes.gzip'.format(dom_dir))
        pd.DataFrame(None, index=cts).to_csv('{}/{}.tsv'.format(dom_dir, dom.lower()), sep='\t')
        
    # whole mouse /human
    dom_dir = '../preset_targets_new/{}'.format(sp)
    os.makedirs(dom_dir, exist_ok=True)

    # write final files for whole mouse and human
    expr = expr.reindex(order, axis=1, level=0)
    expr_h1 = expr.iloc(1)[:int(expr.shape[1]/2)].to_pickle('{}/expression_h1.gzip'.format(dom_dir))
    expr_h2 = expr.iloc(1)[int(expr.shape[1]/2):].to_pickle('{}/expression_h2.gzip'.format(dom_dir))    
    
    markergenes = markergenes.reindex(order, axis=1, level=1)
    markergenes.to_pickle('{}/markergenes.gzip'.format(dom_dir))
    pd.DataFrame(None, index=order).to_csv('{}/{}.tsv'.format(dom_dir, sp), sep='\t')
    pd.concat(colors).to_csv('{}/colors.tsv'.format(dom_dir), sep='\t', header=True)