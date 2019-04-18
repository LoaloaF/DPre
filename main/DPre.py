import pandas as pd
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

import _logger
import _dpre_util as util
import config
from _format_input import TARGET
from targets import Targets
from drivers import Drivers





'''original single inhibitors'''
expr = pd.read_csv('inh/results/genes_cpm_expression.tsv', sep='\t', 
                   index_col='ensg')
expr = expr.loc[expr.index.values[1::2], :]
order = ['BIX', 'Bromospormine', 'CPTH2', 'Chaetocin', 'DZNep', 
         'Garcinol', 'IOX1', 'MS023', 'OICR', 'Sirtinol', 'Tranylcypromine', 
         'UNC0379', 'Valproate', 'Twoi']
expr = expr.reindex(['mean_'+ o for o in order], axis=1)
expr.columns = order
expr.sort_index(axis=1, inplace=True)
# diff = fmt_deseq2_out('inh/significance/deseq2/results/_sig/Twoi_sig')
dirr = 'emb/emb_inh/deseq2/res'
dirr = 'lineages/cell_lineages'
dirr2 = 'lineages/embryo'
d = Drivers(expression=expr, diff_genes='inh/significance/deseq2/results/_sig/Twoi_sig', name='epigenetic inhibitors', ctrl='Twoi', override_diff_names=True)
d.set_colors(config.default_colors)
# d = Drivers(expression=expr, diff_genes=dirr2, name='epigenetic inhibitors', ctrl='Twoi')
# # d.add_expression(expr, new_elements=['Twoi'], check_alignment=False, row_wise_variance=False, show_hist=False)
# d.set_colors(colors=config.default_colors)
order = [
        'DZNep', 'Valproate', 'Bromospormine',  'MS023', 'Tranylcypromine', 'IOX1',  
         'UNC0379',  'CPTH2', 'OICR', 'Sirtinol', 'Garcinol', 'Chaetocin', 'BIX', 
         'Twoi']
d.reorder(order)
# sys.exit( )
# d = d.slice_elements(['DZNep', 'Twoi'])

# T = TARGET('all')
# T.overlap_drivers(d, 'euclid')
# sys.exit()


# t= TARGET('all', False, True)
'''original multi inhibitors'''
expr = pd.read_csv('test_data/genes_ntc_expression.tsv', sep='\t', 
                   index_col='ensg')
order = ['4C-8C DMIT', '4C BMIT', '4C DBMIT', '4C DBMI','4C DBMT',
         '8C DMVIT', '8C DMVI', '8C DMVT', '8C DVIT', 
         'DIB', 'DMB', 'DMI', 'DMVIBTCOUCGSB', 'DMV',  'DVB', 'DVI', 
         'ICM DITCO', 'ICM DTCO', 'ICM ITCO', 
         'MIB', 'MVB', 'MVI', 'VIB', 'none']
expr = expr.reindex(['mean_ESC '+ o for o in order], axis=1)
expr.columns = order
# diff = fmt_deseq2_out('multi_inh/multi_single/deseq2/multi')


c = Drivers(diff_genes='DPre/test_data/deseq2', name='epigenetic inhihbitor combinations',
            expression=expr, ctrl='none', override_diff_names=True)
t = TARGET('embryonic')
t.plot_overlap_distances(c, 'euclid', proportional=False,  

              reorder_to_max_possible_bar=False, 
              cluster_targets=True, 
              cluster_drivers=False,
              
              heatmap_width=1, 
              heatmap_height=1, 
               
              heatmap_range_from_config=False,
              max_possible_bar_range_from_config=False,

              show_max_possible_bar=True, 
              show_target_dendrogram=True,
              show_driver_dendrogram=True, 

              drivers_colorbar = True, 
              targets_colorbar = True,
)
sys.exit()

# c = MeasuredDriverCombinations(diff_genes='multi_inh/multi_single/deseq2/multi', 
#                      expression=expr,  name='epigenetic inhibitor combinations', 
#                      derived_from_drivers=d, composition_mapping=combs, 
#                      ctrl='none', override_diff_names=True)
order = ['4C-8C DMIT', '8C DMVI', '8C DMVT', '8C DVIT', '8C DMVIT',
         'ICM DITCO', 'ICM DTCO', '4C DBMIT', '4C DBMI','4C DBMT',
         'DMV', 'DIB', 'DMB', 'DMI',  'DVB', 'DVI',
         'ICM ITCO', '4C BMIT', 
         'MIB', 'MVB', 'MVI', 'VIB', 
         'none']
c = c.slice_elements(order)


D_3 = dict.fromkeys(['DMV', 'DIB', 'DMB', 'DMI',  'DVB', 'DVI'], config.default_colors[4])
D_4 = dict.fromkeys(['4C-8C DMIT', '8C DMVI', '8C DMVT', '8C DVIT', '8C DMVIT', 
                     'ICM DITCO', 'ICM DTCO', '4C DBMIT', '4C DBMI','4C DBMT'], 
                     config.default_colors[5])
noD_3 = dict.fromkeys(['MIB', 'MVB', 'MVI', 'VIB'], config.default_colors[10])
noD_4 = dict.fromkeys(['ICM ITCO', '4C BMIT'], config.default_colors[11])
print(config.default_colors)
c.set_colors({**D_3, **D_4,  **noD_4,  **noD_3, })


# c._generated = c._make_generated_equivalent()

t = TARGET('embryonic', colors_from_file=True) 
genes = ['Duxf3', 'Gata2', 'Zscan4a', 'Zscan4b', 'Zscan4c', 'Zscan4d', 'Zscan4e', 'Zscan4f', 'Zscan4a', 'Hsp90aa1', 'Hsp90ab1', # 2C genes
    'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d', 'Utf1', 'Dppa2', 'Dppa3', 'Dppa4', 'Il6st', 'Foxo1', 'Prdm4', 'Zic1', 'Il6st',
    'Epcam', 'Tdgf1', 'Cldn3', 'Cldn4', 'Ocln', 'Alpl', 'Crb3', 'Cldn11', 'Inadl','Terf1', # 'Initiation'
    "Gdf3", "Nanog", "Zfp42", "Esrrb", "Prdm1", # 'Maturation'
    "Utf1", "Dppa2", "Dppa3", "Dppa4", "Pecam1`", "Dnmt3l", "Klf2", "Klf5", # 'Stabilization'
    'Nr6a1', 'Tbx21', 'Cldn11', 'Six4', "Prx",
    "Tbx3", 'Nr5a1', 'Nr5a1', 'Nr5a3', 'Nr0b1', 'Gdf3', 'Nr6a1', 'Utf1', 'Lin28a', 'Sall4', 'Fbxo15', 'Snai2',
    'Pkp1', 'Pkp2', 'Lox', 'Serpine1', 'Thbs1', 'Col6a1', 'Tgfb3', 'Tcf4', 'Ocln', 'Prdm1', 'Tgfb1', 'Klf8', 'Zeb1', 'Zeb2', 'Fn1',
    'Pecam1', 'Nr6a1', 'Nr5a1', 'Nr5a2', 'Wnt5a', 'Wnt3a', 'Trp63', 'Trp53', 'Bmi1', 'Stat3', 'Lifr', 'Bcl3', 'Socs3', 'Bmp1ra', 'Cdc20', 'Ccnf',
    'Stat3', 'Fgf4', 'Dnmt3b', 'Notch1', 'Grb2', 'Tcf3', 'Zic2', 'Sall1', 'Lin28a',
    "Ncor1", "Ncor2", "Snai1", "Thy1", "Acta2", "Acta1", 'Actg1', 'Actg2', "Gapdh", "Actb", "Eef1a1",
    "Pou5f1", "Sox2", "Nanog", "Cdh1", "Cdh2", "S100a4", "Esrrb", "Klf4", "Myc", "Lin28a",
    "Pecam1", "Dppa5a", 'Dppa2',
    "Rcor2", "Lin28b", "Foxo4", "Ccnd3", "Ccnb1", "Cbx1", "Dab2", "Gli1", "Tbx2", "Tgfb2",
    ]
genes = pd.Index(genes).drop_duplicates().tolist()

cols_init = dict.fromkeys(['Duxf3', 'Gata2', 'Zscan4a', 'Zscan4b', 'Zscan4c', 'Zscan4d', 'Zscan4e', 'Zscan4f', 'Zscan4a', 'Hsp90aa1', 'Hsp90ab1', # 2C genes
    'Mef2a', 'Mef2b', 'Mef2c', 'Mef2d', 'Esrrb', 'Utf1', 'Dppa2', 'Dppa3', 'Dppa4', 'Il6st', 'Foxo1', 'Prdm4', 'Zic1', 'Il6st',
    'Epcam', 'Tdgf1', 'Cldn3', 'Cldn4', 'Ocln', 'Alpl', 'Crb3', 'Cldn11', 'Inadl','Terf1'], config.default_colors[15])
cols_mat = dict.fromkeys(["Gdf3", "Nanog", "Zfp42", "Esrrb", "Prdm1"], config.default_colors[16])
cols_stab = dict.fromkeys(["Klf2", "Klf5", 'Nr6a1', 'Tbx21', 'Cldn11', 'Six4', "Prx",
    "Tbx3", 'Nr5a1', 'Nr5a1', 'Nr5a3', 'Nr0b1', 'Gdf3', 'Nr6a1', 'Utf1', 'Lin28a', 'Sall4', 'Fbxo15', 'Snai2',
    'Pkp1', 'Pkp2', 'Lox', 'Serpine1', 'Thbs1', 'Col6a1', 'Tgfb3', 'Tcf4', 'Ocln', 'Prdm1', 'Tgfb1', 'Klf8', 'Zeb1', 'Zeb2', 'Fn1',
    'Pecam1', 'Nr6a1', 'Nr5a1', 'Nr5a2', 'Wnt5a', 'Wnt3a', 'Trp63', 'Trp53', 'Bmi1', 'Stat3', 'Lifr', 'Bcl3', 'Socs3', 'Bmp1ra', 'Cdc20', 'Ccnf',
    'Stat3', 'Fgf4', 'Dnmt3b', 'Dnmt3l', 'Notch1', 'Grb2', 'Tcf3', 'Zic2', 'Sall1', 'Lin28a',
    "Ncor1", "Ncor2", "Snai1", "Thy1", "Acta2", "Acta1", 'Actg1', 'Actg2', "Gapdh", "Actb", "Eef1a1",
    "Pou5f1", "Sox2", "Nanog", "Cdh1", "Cdh2", "S100a4", "Esrrb", "Klf4", "Myc", "Lin28a",
    "Pecam1", "Dppa5a", 'Dppa2', 
    "Rcor2", "Lin28b", "Foxo4", "Ccnd3", "Ccnb1", "Cbx1", "Dab2", "Gli1", "Tbx2", "Tgfb2"], config.default_colors[17])

t.map_overlap(c, 'euclid', proportional=False, differential=True, 
              filename='combinations_2C.geneset_all.pdf',

              gene_number=None, 
              custom_target_genelist=genes,
              specific_markergenes = None,
                specific_markergenes_show_none=True,
              
              reorder_to_max_possible_bar=False, 
              cluster_genes=True, 
              cluster_drivers=False,
              
              heatmap_width=1, 
              heatmap_height=1, 
               
              show_genelabels=True,  
              genelabels_size=.9, 

              sum_plot_central_metric='mean',
              heatmap_range_from_config=True,
              sum_plot_xlim_from_config=True, 
              max_possible_bar_range_from_config=True,

              show_max_possible_bar=True, 
              show_sum_plot=True,
              show_gene_dendrogram=True,
              show_driver_dendrogram=True, 

              drivers_colorbar = True, 
            #   genes_colorbar = {**cols_stab, **cols_init, **cols_mat}, 
              genes_colorbar = False,


)

# test = c.slice_elements(['8C DMVI', '4C DBMI'], 'test-combination')
# c_threes = c.slice_elements(['DMV', 'DMI', 'DVI', 'MVI', 'DMB',  'DVB', 'DIB', 'MVB', 
#                              'MIB', 'VIB'], '3 inhi.')
# test._generated = test._make_generated_equivalent()

# c_threes._generated = c_threes._make_generated_equivalent()
# t = TARGET('embryonic', colors_from_file=True) 
# c._get_explained_effect('euclid', trg=t)
# t = TARGET('embryonic', colors_from_file=True) 
# c._get_explained_effect('euclid', trg=t)
# t = TARGET('mesoderm', colors_from_file=True) 
# c._get_explained_effect('euclid', trg=t)


# valid = ['all', 'blood mesoderm', 'embryonic', 'endoderm', 'germ cells',
#          'mesoderm', 'neural crest', 'neuroectoderm', 'surface ectoderm']
# genes = c._detec_genes.intersection(c._generated._detec_genes)
# for v in valid:

#     t = TARGET(v)
#     match = len(t._diff_mgs.index.intersection(genes))
#     ofall = len(t._diff_mgs.index)
#     print('{}: {}/{}, {}%'.format(v, match, ofall, match/ofall))






# t = TARGET('all', False, True)
# sl = ['B cell (active)', 'B cell (germinal center)', 'B cell (naive)', 'B cell (resting)', 
# 'B cells CD19+', 'B cells CD43-', 'BMDM', 'BMDM +CpolyG', 'BMDM +IFNa', 'BMDM +IFNg', 'BMDM +IL-1b', 'BMDM +IL-4', 'BMDM +MALP2', 'BMDM +TGFb', 'BMDM +lipidA 360 mins', 'BMDM -lipidA', 'Basophilic erythroblast', 'Bone marrow', 'Bone marrow macrophages', 'CD4+ T NKG2D+ cells', 'CD4+ T NKG2D- cells', 'CD4+ T cells', 'CD4+ T cells +IL-21', 'CD4+ T cells -IL-21', 'CD4+ Th1 cells', 'CD4+ Th17 cells ', 'CD4+ Th2 cells', 'CD4+ Th9 cells', 'CD4+ Treg cells', 'CD4+ iTreg cells', 'CD4T DP CD3hi', 'CD4T DP CD3lo', 'CD8+ T cells', 'CD8+ T cells +IL-2', 'Common lymphoid progenitor', 'Common myeloid progenitor', 'DC CD11+', 'DC CD8+', 'DN1e', 'Dendritic cells (BMDC)', 'Double negative 1 (DN1) thymic cells', 'Double negative 2a (DN2a) thymic cells', 'Double negative 2b (DN2b) thymic cells', 'Double negative 3 (DN3) thymic cells', 'Double positive thymic cells', 'Eosinophils -IL-10 +LPS', 'Eosinophils -IL-10 -LPS', 'Erythroblast', 'Erythrocytes type A', 'Erythrocytes type B', 'Erythroid cells (TER119+)', 'Erythroid progenitor', 'Granulocyte-monocyte progenitor', 'Granulocytes', 'Hematopoietic stem cells', 'Hematopoietic stem cells (MPP1)', 'Hematopoietic stem cells (MPP2)', 'Hematopoietic stem cells (MPP3)', 'Hematopoietic stem cells (MPP4)', 'Innate lymphoid cell (Type 2)', 'Invariant NKT', 'Long-term HSC', 'Macrophages +IL10 +LPS', 'Macrophages +IL10 -LPS', 'Macrophages -IL10 +LPS', 'Macrophages -IL10 -LPS', 'Mast cell +IL10 +LPS', 'Mast cell +IL10 -LPS', 'Mast cell -IL10 +LPS', 'Mast cell -IL10 -LPS', 'Mast cells (CD25+)', 'Mast cells (CD25-)', 'Megakaryoblast', 'Megakaryocyte erythroid progenitor', 'Microglia', 'Monocyte', 'Monocytes', 'Multipotent progenitor', 'Myelocytes', 'Myeloid cells', 'NKT (derived) ', 'NKT stage1', 'NKT stage2', 'NKT stage3', 'Naive CD4+ T cells, anti-CD3-28 treated', 'Natural killer cells', 'Neutrophil +IL10 +LPS', 'Neutrophil +IL10 -LPS', 'Neutrophil -IL10 +LPS', 'Neutrophil -IL10 -LPS', 'Orthochromatic erythroblast', 'Polychromatic erythroblast', 'Proerythroblast', 'Promyelocytes', 'Spleen', 'Thymus', 'pre-B cells', 'pro-B (fraction B)', 'pro-B (fraction CC`)', 'pro-B cells', 'splenic B cells', 'splenic DC +IL10 +LPS', 'splenic DC +IL10 -LPS', 'splenic DC -IL10 +LPS', 'splenic DC -IL10 -LPS', 'Caecum (E13.5)', 'Colon', 'Colon epithelial cells', 'Ileum', 'Intestinal stem cell', 'Intestine', 'Large intestine', 'Liver', 'Lung', 'Medullary thymic epithelial', 'Medullary thymic epithelial cell (immature)', 'Medullary thymic epithelial cell (mature)', 'Pancreatic cell (alpha)', 'Pancreatic cell (beta)', 'Prostate basal cells', 'Prostate luminal cells', 'Small intestine cells', 'Mature sperm', 'Spermatids', 'Spermatids (elongated)', 'Spermatids (round)', 'Spermatocytes', 'Spermatogonia', 'Testicular cells', 'Adipocytes (brown)', 'Adipocytes (cultured, brown)', 'Adipocytes (cultured, white)', 'Adipocytes (white)', 'Adipose', 'Adrenal medulla', 'Atrioventricular canal (E12.5)', 'Carotid body', 'Chondrocyte (Rib)', 'Cortical thymic epithelial', 'Dermal fibroblast', 'Embryonic kidney fibroblast', 'Epididymis', 'Fatpad', 'Heart', 'Intestinal subepithelial myofibroblasts (adult)', 'Intestinal subepithelial myofibroblasts (neonatal)', 'Kidney', 'Lung fibroblasts', 'MEF', 'MSC (bone)', 'MSC (skin)', 'MSC (thymus)', 'Mammary tissue', 'Myoblast', 'Neonatal tail fibroblast', 'Oviduct', 'Skeletal muscle', 'Somatic gonadal cell', 'Striated muscle', 'Tail fibroblast', 'Uterus', 'Vas deferens', 'pre-Adipocytes (brown)', 'pre-Adipocytes (white)', 'Caudal brain neural epithelium (E8.5)', 'Central neural epithelium (E10.5)', 'Cranial mesenchyme (E9.5)', 'Dorsal control neural epithelium (E10.5)', 'Epidermal ectoderm (E9.5)', 'Floor plate (E8.5)', 'Lateral nasal prominence (E10.5)', 'Lateral prominence neural epithelium (E10.5)', 'Lower molar', 'Mandibular arch (E10.5)', 'Mandibular arch (E9.5)', 'Maxillary arch (E10.5)', 'Maxillary arch (E9.5)', 'Mesoderm (E8.5)', 'Non-floor plate neural epithelium (E8.5)', 'Olfactory pit (E10.5)', 'Upper molar', 'Astrocyte', 'Brain', 'Cerebellar granular neurons', 'Cerebellum', 'Cortex', 'Cortex neuron', 'Cortical neurons (E16.5)', 'Dentate gyrus', 'Dorsal root ganglia', 'Eye', 'Frontal cortex', 'Hippocampus', 'Hippocampus neuron', 'Hippocampus tissue', 'Medial ganglionic eminence neurons', 'Motor neurons', 'Neocortex cortical plate', 'Neocortex subventricular zone', 'Neocortex ventricular zone', 'Neural progenitor cell', 'Neuron', 'Nucleus accumbens', 'Olfactory bulb', 'Oligodendrocyte (myelinating)', 'Oligodendrocyte (newly formed)', 'Oligodendrocyte (precursor)', 'Pineal gland', 'Retina', 'Spinal cord', 'Striatum', 'Substantia nigra', 'Subventricular zone', 'Telencephalon', 'Ventral tegmental', 'White-matter glia', 'Zona limitans intrathalamica', 'Basal epidermis (E14.5)', 'Dermis (E18.5)', 'Epidermis (E18.5)', 'Hair follicle', 'Hair follicle stem cells', 'Hair follicle transit amplifying cells', 'Keratinocyte', 'Skin (E16.5)', 'Skin epithelial', 'Suprabasal epidermis (E14.5)']
# t = t.slice_elements(sl)
# t = TARGET('embryonic', False, True)

# t.plot_prediction_accuracy(c, 'euclid', cluster_targets=False, targets_colorbar=True, heatmap_range_from_config=True,
# max_possible_bar_range_from_config=False, show_targetlabels=True, reorder_to_max_possible_bar='combination', filename='ctrl difference sorted acc.png')






# t._get_from_overlap(d, 'euclid', True)

# sys.exit()
# agg_mean.unstack((0,1)).reindex(agg_mean.index.unique(2))
# c.add_expression(expr, show_hist=False, check_alignment=False)
# col = dict.fromkeys(['8C_DMVIT', '4C-8C_DMIT', '8C_DMVT', '8C_DVIT', '8C_DMVI',
#                     '4C_DBMIT', '4C_DBMT', '4C_DBMI', 'ICM_DITCO', 'ICM_DTCO', 'DMV',
#                     'DMI', 'DMB', 'DVI', 'DVB', 'DIB', 'DMVIBTCOUCGSB'], '#00806c')
# col2 = dict.fromkeys(['4C_DBMIT', '4C_DBMT', '4C_DBMI', 'DMB', 'DVB', 'DIB', 
#                       'DMVIBTCOUCGSB', '4C_BMIT', 'MIB', 'MVB', 'VIB'], '#5a3500')
# col3 = dict.fromkeys( ['4C_BMIT', 'MIB', 'MVB', 'VIB'], '#d90083')
# col.update(col2)
# col.update(col3)
# c.set_colors(col)
# # order = ['DMV', '8C_DMVI', 'none_ctrl']
# # order = ['DMV', 'none_ctrl'] 
# order = ['DMV', '8C_DMVIT', '8C_DMVI', '4C-8C_DMIT', '8C_DMVT', 'ICM_DTCO', '8C_DVIT', 'ICM_DITCO', 
#          '4C_DBMI', '4C_DBMIT', '4C_DBMT', 'DVB', 'DMB', 'DIB', 'DMI', 
#          'DVI', 'MVB', 'MVI', 'MIB', 'VIB', '4C_BMIT', 'ICM_ITCO', 'none_ctrl']
# c = c.slice_elements(order)
# c.convert_to_drv_diff()
# # print(c.names)
# # c.reorder(order)


# # c.check_prediction()
# # sys.exit()



# '''single inhibitors'''
# # expr = pd.read_csv('multi_inh/multi_single/rsem-genes/genes_ntc_expression.tsv', sep='\t', 
# #                    index_col='ensg')
# # order = ['BIX', 'Bromospormine', 'CPTH2', 'Chaetocin', 'DZNep', 
# #          'Garcinol', 'IOX1', 'MS023', 'OICR', 'Sirtinol', 'Tranylcypromine', 
# #          'UNC0379', 'Valproate', 'Twoi']
# # expr = expr.reindex(['mean_ESC_'+ o for o in order], axis=1)
# # diff = fmt_deseq2_out('multi_inh/multi_single/deseq2/single')
# # d = drivers(diff, down_markergenes=False, name='sepigenetic inhibitors', ctrl='Twoi')
# # d.add_expression(expr, new_elements=['Twoi'], show_hist=False)
# # order = ['BIX', 'CPTH2', 'Chaetocin',  
# #          'Garcinol', 'IOX1',  'OICR', 'Sirtinol', 'Tranylcypromine', 
# #          'UNC0379', 'DZNep', 'Valproate', 'MS023', 'Bromospormine', 'Twoi']
# # d.reorder(order)
# # d = d.slice_elements(['DZNep', 'MS023', 'Valproate', 'Twoi'])





# '''multi inhibitors'''

# # expr = pd.read_csv('multi_inh/multi_single/rsem-genes/genes_ntc_expression.tsv', sep='\t', 
# #                    index_col='ensg')
# # order = [ 'none', '4C-8C_DMIT', '4C_BMIT', '4C_DBMI', '4C_DBMIT', '4C_DBMT',
# #          '8C_DMVI', '8C_DMVIT', '8C_DMVT', '8C_DVIT', 'ICM_DITCO', 'ICM_DTCO', 
# #          'ICM_ITCO', 
# #          'DIB', 'DMB', 'DMI', 'DMV', 'DMVIBTCOUCGSB', 'DVB', 'DVI', 'MIB', 
# #          'MVB', 'MVI', 'VIB']
# # expr = expr.reindex(['mean_ESC_'+ o for o in order], axis=1)
# # expr = expr.reindex(columns=expr.columns.sort_values())
# # diff = fmt_deseq2_out('multi_inh/multi_single/deseq2/multi')
# # d = drivers(diff, down_markergenes=True, name='epigenetic inhibitor combinationshm', ctrl='none_ctrl')
# # d.add_expression(expr, show_hist=False)
# # col = dict.fromkeys(['8C_DMVIT', '4C-8C_DMIT', '8C_DMVT', '8C_DVIT', '8C_DMVI',
# #                     '4C_DBMIT', '4C_DBMT', '4C_DBMI', 'ICM_DITCO', 'DMV',
# #                     'DMI', 'DMB', 'DVI', 'DVB', 'DIB', 'DMVIBTCOUCGSB'], '#00806c')
# # col2 = dict.fromkeys(['4C_DBMIT', '4C_DBMT', '4C_DBMI', 'DMB', 'DVB', 'DIB', 
# #                       'DMVIBTCOUCGSB', '4C_BMIT', 'MIB', 'MVB', 'VIB'], '#5a3500')
# # col.update(col2)
# # d.set_colors(col)
# # order = ['DMV', 'none_ctrl']
# # d = d.slice_elements(order)




# # d = np.log2(expr.loc(0)['ENSMUSG00000057666'])
# # plt.plot(d)
# # sys.exit(0)
    
# # t.map_overlap(d, 'euclid', difference_to='none_ctrl', proportional=True, set_not_detected_to_zero=False)
# # t.map_overlap(d, 'euclid', difference_to='Twoi', cluster_columns=True, show_max_possible_bar=False, proportional=True, colorbar_genes=None,
# #     show_left_dendogram=True, show_top_dendogram=False, set_not_detected_to_zero=False,
# #     show_not_intersecting_genes=True, show_sum_plot=False)

'''embryo weird'''
# expr = pd.read_csv('emb/emb_inh/rsem_res/genes_cpm_expression_cut2.tsv', sep='\t', 
#                    index_col='ensg')
# expr = expr.iloc[1::2]
# # order = ['BIX', 'Bromospormine', 'CPTH2', 'Chaetocin'l/;, 'DZNep', 
# #          'Garcinol', 'IOX1', 'MS023', 'OICR', 'Sirtinol', 'Tranylcypromine', 
# #          'UNC0379', 'Valproate', 'Twoi']
# [expr.drop(c, axis=1, inplace=True) for c in expr.columns if 'mean_' not in c]
# expr = expr.reindex([
#        'mean_embryo 2C',
#        'mean_embryo 4C',
#        'mean_embryo 8C',
#       'mean_embryo ICM',
#             'mean_Twoi',
#   'mean_embryo early2C',
# 'mean_embryo miioocyte',
#    'mean_embryo zygote'], 
#    axis=1)
# expr.columns = ['2C', '4C', '8C',  'ICM', 'Twoi', 'early2C', 'miioocyte', 'zygote']
# dirr = 'emb/emb_inh/deseq2/res'
# t = Targets(diff_genes=dirr, expression=expr, name=None, ctrl='Twoi', override_diff_names=True, use_down_mgs=True)
# t.reorder(['miioocyte', 'zygote', 'early_2C', '2C', '4C', '8C', 'ICM'])
# t.set_colors(config.default_colors)

# t = t.slice_elements(['miioocyte', 'zygote', 'early_2C', '2C', '8C', 'ICM'])
# t = t.slice_elements(['miioocyte', 'zygote'])
# t.add_expression(expr, new_elements=['Twoi'], row_wise_variance=False, check_alignment=False)
# t = t.slice_elements([ 'ICM', 'early_2C',   'miioocyte', 'zygote', '2C', '4C', '8C'])





# t.overlap_drivers_data(d.expr.xs('z', axis=1, level=1), 'euclid',
#                        d.not_expr_zvalue)

# o = t.get_from_overlap(d, 'intersect', pd.IndexSlice[:, 'ICM', :], differential=True, standardize_down_mgs=True, eff_proportional=False)

# print('\n\n')
# print(o[0])
# print(o[1])
# print(o[2])
# print(o[3])..


# t = TARGET('blood mesoderm', False)

# os.chdir('DPre/default_targets')
# for f in glob.glob('*.tsv'):
#     names = pd.read_csv(f, sep='\t', header=None).iloc(1)[1].values
#     get = t.expr.loc(1)[(names, 'z')]
#     get.columns = get.columns.unique(0)
#     get.to_csv(f[:-4]+'/z_expression.tsv', sep='\t')
    # z_expr = pd.read_csv(f[:-4]+'/z_expression.tsv',
    #                      sep = '\t',S
    #                      index_col='ensg')
                         


# t.plot_mapped_expression(t)
# t = t.slice_elements(['4C'], log=False)


# t = TARGET('all', True)
# # d = t.plot_overlap_distances(d, 'euclid', differential=True, title=False, show_targetlabels=True, proportional=False, show_driverlabels=True,
#                              drivers_colorbar=False, cluster_targets=True, show_driver_dendrogram=True, distance_bar_range_from_config=True,
#                              cluster_drivers=True, reorder_to_max_possible_bar=False, show_max_possible_bar=False, show_colorbar_legend=True,
#                              heatmap_range_from_config=False, targets_colorbar=True, heatmap_width=1,
#                              filename='normed_all_multinba.png')[2]
# print(d['up'][0].columns.tolist())

# t = t.slice_elements(['2C Embryo', '4C Embryo'])


# t.overlap_drivers(d, 'euclid')






# t = t.slice_elements(['2C Embryo'])
# t.plot_overlap_distances(d, 'euclid', differential=True, title=False, show_targetlabels=True, proportional=False, show_driverlabels=True,
#                              drivers_colorbar=False, cluster_targets=True, show_driver_dendrogram=False, distance_bar_range_from_config=True,
#                              cluster_drivers=False, reorder_to_max_possible_bar=False, show_max_possible_bar=True, show_colorbar_legend=True,
#                              heatmap_range_from_config=True, targets_colorbar=True, heatmap_width=1,  
#                              filename='multi_eucl_emb.png')[2]
# order = ['Neonatal tail fibroblast', 'Dermal fibroblast', 'Tail fibroblast', 'Lung fibroblasts', 'Myoblast', 'MEF', 'pre-Adipocytes (white)', 
#          'pre-Adipocytes (brown)', 'MSC (skin)', 'MSC (bone)', 'MSC (thymus)', 'Adipocytes (cultured, brown)', 'Chondrocyte (Rib)',
#          'Adipocytes (cultured, white)', 'Intestinal subepithelial myofibroblasts (adult)', 'Intestinal subepithelial myofibroblasts (neonatal)']
# t = TARGET('all', True)
# t = t.slice_elements(order) 
# t.plot_overlap_distances(d, 'euclid', differential=True, title=False, show_targetlabels=True, proportional=False, show_driverlabels=False,
#                              drivers_colorbar=False, cluster_targets=True, show_driver_dendrogram=True, distance_bar_range_from_config=True,
#                              cluster_drivers=False, reorder_to_max_possible_bar=False, show_max_possible_bar=True, show_colorbar_legend=False,
#                              heatmap_range_from_config=True, targets_colorbar=True, heatmap_width=1,
#                              filename='multi_eucl_meso.png')[2]
# t = TARGET('blood mesoderm', True) 
# order = ['BMDM -lipidA', 'BMDM +TGFb', 'BMDM', 'BMDM +IL-4', 'BMDM +lipidA 360 mins', 'BMDM +CpolyG', 'BMDM +IFNa', 'Macrophages +IL10 -LPS',
#          'Macrophages -IL10 -LPS', 'BMDM +IL-1b', 'BMDM +MALP2', 'Dendritic cells (BMDC)', 'Myeloid cells', 'BMDM +IFNg', 
#          'Macrophages +IL10 +LPS', 'Macrophages -IL10 +LPS', 'Mast cell +IL10 +LPS', 'Mast cell -IL10 +LPS', 
#          'Eosinophils -IL-10 +LPS', 'Eosinophils -IL-10 -LPS', 'Myelocytes', 'DC CD11+', 'Neutrophil +IL10 -LPS', 'Neutrophil -IL10 -LPS',
#          'Neutrophil +IL10 +LPS', 'Neutrophil -IL10 +LPS', 'Monocyte', 'pro-B cells', 'Hematopoietic stem cells (MPP3)', 
#          'Hematopoietic stem cells (MPP4)', 'Hematopoietic stem cells', 'Mast cells (CD25+)', 'Mast cells (CD25-)', 'Promyelocytes', 
#          'Bone marrow macrophages', 'Megakaryoblast', 'Granulocytes', 'splenic DC +IL10 +LPS', 'splenic DC -IL10 +LPS', 
#          'Mast cell +IL10 -LPS', 'Mast cell -IL10 -LPS', 'Monocytes', 'DC CD8+', 'splenic DC +IL10 -LPS', 'splenic DC -IL10 -LPS']
# order = ['Double negative 2a (DN2a) thymic cells', 'Double negative 2b (DN2b) thymic cells', 
#          'Double negative 3 (DN3) thymic cells', 'Double positive thymic cells', 'CD4+ Th17 cells ', 
#          'CD4+ Th9 cells', 'DN1e', 'CD4T DP CD3hi', 'CD4T DP CD3lo', 'CD4+ Treg cells', 
#          'Naive CD4+ T cells, anti-CD3-28 treated', 'NKT stage3', 'NKT stage1', 
#          'NKT stage2', 'CD4+ T NKG2D+ cells', 'CD4+ T NKG2D- cells', 
#          'CD4+ T cells +IL-21', 'CD4+ T cells -IL-21', 'CD4+ Th1 cells', 'CD4+ Th2 cells']
# t = t.slice_elements(order)
# d.plot_mapped_expression(t, sort=True, colored=True)
 
# t = TARGET('embryonic', True) 
# t.plot_prediction_accuracy(c, 'euclid', title=False, show_targetlabels=True, proportional=False, show_driverlabels=True,
#                              drivers_colorbar=False, cluster_targets='accuracy', max_possible_bar_range_from_config=True, 
#                              heatmap_range_from_config=False, targets_colorbar=True, heatmap_width=1, show_single_drivers=True,
#                              filename='multi_eucl_blmeso.png')
# t.plot_overlap_distances(d, 'euclid', differential=True, title=False, show_targetlabels=True, proportional=False, show_driverlabels=True,
#                              drivers_colorbar=False, cluster_targets=True, show_driver_dendrogram=False, distance_bar_range_from_config=False,
#                              cluster_drivers=False, reorder_to_max_possible_bar=False, show_max_possible_bar=True, show_colorbar_legend=True,
#                              heatmap_range_from_config=False, targets_colorbar=True, heatmap_width=1,


# t.map_overlap(d, 'euclid', proportional=False, gene_number=30, show_genenames=True,  drivers_colorbar=True, differential=True, cluster_genes=True, genenames_size=1.3,
#               show_max_possible_bar=False, reorder_to_max_possible_bar=True, heatmap_width=2, specific_genes_only=True, show_sum_plot=True, filename='yead2.pdf', distance_bar_range_from_config=True,
#               cluster_drivers=True, specific_genes=['Dnmt3l', 'Dppa3', 'Dnmt3b', 'Satb2',
#                     'Csnk2a1', 'Zfp57', 'Prkaa2', 'Aurkc', 'Pcgf1', 'Ctcfl', 'Sfmbt2', 'Taf9', 'Brdt'])
            #   genes_colorbar= {**dict.fromkeys(['Dynap', 'Gm2016', 'Bhmt'], '#ff7f00'), **{'Rbm44': '#a65628'} })
# pd.DataFrame().apply


# t.driver_prediction(d, 'euclid', colored=True, load_prediction=False, off_target_threshhold=.7, colorbar=True, show_max_only=True, 
# #                     exclude_single_predictions=True, flip=False, plot_celltype_identity=True, filename='dinal.png')
# t.driver_prediction(d, 'euclid', driver_number_threshold=None, driver_effect_threshold=None, alt_comb_numbers=20, show_driver_colorlegend=True,
#                     show_max_only=True, off_target_numbers=2, proportional=False, alt_combs_threshhold=None, show_singlebars=True,
#                     exclude_single_predictions=False, flip=False, load_prediction=False, convert_to_proportional=False, show_driverlabels=False,
#                     colorbar=True, colored_drivers=True, show_updown=False, filename='DMV_all_prop2.png', colored_targets=False,
#                     off_target_threshhold=None, specify_comb=['DZNep', 'Valproate', 'MS023'], plot_celltype_identity=True)

# t.plot_driver_effects(d, 'euclid',  colorbar=False, show_target_labels=True, show_updown=True, 
#                 target_numbers=50, colored=True, plot_celltype_identity=False, target_threshhold=None, 
#                 proportional=False, driver_numbers=3, plot_layout=(1,3), xlim_from_config=True,
#                 filename='testop4.png') 


'''multi inhibitors'''
# d.add_expression(expr, False, False)
# order = ['Twoi_ctrl', 'BIX', 'Bromospormine', 'CPTH2', 'Chaetocin', 'DZNep',  
#          'Garcinol', 'IOX1', 'MS023', 'OICR', 'Sirtinol', 'Tranylcypromine', 
#          'UNC0379', 'Valproate', 
#          'none_ctrl', '4C-8C_DMIT', '4C_BMIT', '4C_DBMI', '4C_DBMIT', '4C_DBMT',
#          '8C_DMVI', '8C_DMVIT', '8C_DMVT', '8C_DVIT', 'ICM_DITCO', 'ICM_DTCO', 
#          'ICM_ITCO', 
#          'DIB', 'DMB', 'DMI', 'DMV', 'DMVIBTCOUCGSB', 'DVB', 'DVI', 'MIB', 
#          'MVB', 'MVI', 'VIB']
# d.reorder(order)
# col_single =  ['#4daf4a'] *13     
# col_fates =  ['#984ea3'] *12
# col_trials =  ['#ff7f00'] *11
# cols = ['#999999'] + col_single + ['#999999'] + col_fates + col_trials
# d.set_colors(cols)




# d.diff_gene_numbers(figsize=(7,4), tight_layout=True, colored=True)
# t.map_expression(d, colored=True, limit_elements=20, sort=True, figsize=(7,5))
# d.plot_mapped_expression(t, sort=True, limit_elements=15, colored=True, figsize=(10,6))



'''embryonic stages'''
# expr = pd.read_csv('emb/original/rsem-genes2/genes_cpm_expression.tsv', sep='\t', index_col='ensg')
# expr.drop(columns=['name', 'strand', 'loc'], inplace=True)
# expr = expr[expr.columns.values[::2]]
# expr = expr.reindex(expr.columns.sort_values(), axis=1)
# diff = fmt_deseq2_out('emb/emb_inh/deseq2/res', p_th=.01)
# t = targets(diff, name='embryonic stages', down_markergenes=False)

# t.add_expression(expr, True)
# t.reorder(['miioocyte', 'zygote', 'early_2C', '2C', '4C', '8C', 'ICM'])

# t.plot_diff_gene_numbers(sort=True, colored=False, tight_layout=True, figsize=(6,4))

# t.plot_overlap_distances(d, which='euclid', colorbar_rows=False, colorbar_columns=False, norm_to='Twoi')







# t.plot_diff_intersection(d)

# t.plot_overlap_distances(d, norm_to='none_ctrl', colored=True, cluster_columns=False, figsize=(20,6))
# t.plot_overlap_distances(d, which='intersect', norm_to='Twoi', colorbar_rows=True, 
#                          show_xlabels=True, cluster_columns=True,
                         
#                          )
# t.plot_overlap_distances(d)
# print(t.expr)
# print(d.expr)

'rapid transcriptional similarity visualizaiton tool for cross datasets predicitons'








'''
z and log2 histogram

        fig, axes = plt.subplots(nrows=2, figsize=(4, 3))
        fig.subplots_adjust(hspace=.5)
        ax = axes
        for i, ax in enumerate(axes):
            if not i:
                m_eff = self._expr.xs('z', 1, 1)
                g_eff = self._generated._expr.xs('z', 1, 1)
                s_eff = self._drivers._expr.xs('z', 1, 1)

                title = ('Distribution of log2- and z-expression between measured combinations, \n'
                         'generated combinations and single inhibitors\n\n\nz-expression')
                ax.set_title(title, pad=5)                 

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlabel('z expression')
                ax.set_ylabel('density')
                # ax.set_ylim(0, 4.5)

                ax.hist(m_eff.values.flatten(), density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled',  label='measured combinations')
                ax.hist(g_eff.values.flatten(), density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled',  label='generated combinations')
                ax.hist(s_eff.values.flatten(), density=True, bins=200, alpha=.6, color='#fc804b', histtype='step',  label='single inhibitors')
                ax.legend(loc='best')
            if i:
                m_eff = self._expr.xs('log2', 1, 1)
                # g_eff = self._generated._expr.xs('log2', 1, 1)
                s_eff = self._drivers._expr.xs('log2', 1, 1)

                title = ('log2-expression')
                ax.set_title(title)                 

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlabel('log2 expression')
                ax.set_ylabel('density')
                # ax.set_ylim(0, 4.5)
                ax.hist(m_eff.values.flatten(), density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled',  label='measured combinations')
                # ax.hist(g_eff.values.flatten(), density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled',  label='generated combinations')
                ax.hist(s_eff.values.flatten(), density=True, bins=200, alpha=.6, color='#fc804b', histtype='stepfilled',  label='single inhibitors')
                plt.legend(loc='best')
        util.save_file(fig, 'log2 z distribution comparison.png')
        # plt.show()




discreapny histogram


        m_ef, g_ef, _, expl, n_exp = data
        m_eff = m_ef.mask(m_ef>2, 2).mask(m_ef<-2, -2)
        g_eff = g_ef.mask(g_ef>2, 2).mask(g_ef<-2, -2)
        expl = expl.mask(expl>2, 2)
        n_expl = n_exp.mask(n_exp>2, 2)

        fig, axes = plt.subplots(nrows=2, figsize=(4, 3))
        fig.subplots_adjust(hspace=.5)

        for i, ax in enumerate(axes):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylim((0,6.5))
            ax.set_xlabel('absolute change in log2-, z-transf. expression')
            ax.set_ylabel('density')
            ax.annotate('>', (.93, -.06), xycoords='axes fraction')

            m_label = 'measured combinations'
            g_label = 'generated combinations'
            d_label = 'discrepancy'
            e_label = 'explained effect'

            if not i:
                title = ('Comparison of `mean` and `sum` single inhibitor combining methods:\n'
                         'discrepancy between measured /comput. generated and explained measured effect\n'
                         '\n\ncombined single inhibitors effects using: `mean`')

                m_data = m_eff.abs().values.flatten()
                g_data = g_eff.abs().values.flatten()
                d_data = n_expl.values.flatten()
                e_data = expl.values.flatten()
                
                ax.annotate('explained max: 72.7', xy=(.031, 6.1), xycoords='data')
                ax.annotate('generated max: 16.8', xy=(.031, 5.8), xycoords='data')
                ax.annotate('measured max: 12.5', xy=(.031, 5.5), xycoords='data')
            elif i:
                ax.annotate('explained max: 65.5', xy=(.031, 6.1), xycoords='data')
                ax.annotate('generated max: 21.2', xy=(.031, 5.8), xycoords='data')
                ax.annotate('measured max: 12.5', xy=(.031, 5.5), xycoords='data')
                title = 'combined single inhibitors effects using: `sum`'

                g_ef = self._make_generated_equivalent(method='sum')._expr_eff
                g_ef = g_ef.reindex(genes).drop('Twoi', axis=1, level=1).xs('z', 1, 2)['up']
                g_ef.fillna(0, inplace=True)
                n_exp = abs(g_ef - m_ef.values)
                exp = (abs(m_ef.values) - n_exp)
                expl = exp.mask(exp < 0, 0)
                
                n_expl = n_exp.mask(n_exp>2, 2)
                expl = expl.mask(expl>2, 2)
                g_eff = g_ef.mask(g_ef>2, 2).mask(g_ef<-2, -2)
                
                m_data = m_eff.abs().values.flatten()
                g_data = g_eff.abs().values.flatten()
                d_data = n_expl.values.flatten()
                e_data = expl.values.flatten()

            ax.set_title(title)

            m_mean = m_data.mean().round(3)
            g_mean = g_data.mean().round(3)
            d_mean = d_data.mean().round(3)
            e_mean = e_data.mean().round(3)

            ax.vlines([m_mean, g_mean, d_mean, e_mean], ymin=0, ymax=4, 
                      color=['purple', '#ffe602', '#ad3232', '#3caf36'], alpha=.6, linewidth=[1.4, 1.4, .4, .4] )
            ax.annotate('mean '+str(m_mean), (m_mean, 3.8), xycoords='data')
            ax.annotate(str(g_mean), (g_mean, 3.2), xycoords='data')
            ax.annotate(str(d_mean), (d_mean, 3), xycoords='data')
            ax.annotate(str(e_mean), (e_mean, 3.8), xycoords='data')
            

            ax.hist(m_data, density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled', range=(0 , 2), label=m_label)
            ax.hist(g_data, density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled', range=(0, 2), label=g_label)
            ax.hist(d_data, density=True, bins=200, alpha=.9, color='#ad3232', histtype='step', range=(0, 2), label=d_label)
            ax.hist(e_data, density=True, bins=200, alpha=.9, color='#3caf36', histtype='step', range=(0, 2), label=e_label)
            if not i:
                ax.legend(loc='best')

        util.save_file(fig, 'combination discrepancy distributions.png', None)
        plt.show()





diff effect histogram
       
        m_ef, g_ef, _, expl, n_exp = data
        m_eff = m_ef.mask(m_ef>2, 2).mask(m_ef<-2, -2)
        g_eff = g_ef.mask(g_ef>2, 2).mask(g_ef<-2, -2)
        expl = expl.mask(expl>2, 2)
        n_expl = n_exp.mask(n_exp>2, 2)

        fig, axes = plt.subplots(nrows=2, figsize=(4, 3))
        fig.subplots_adjust(hspace=.5)

        for i, ax in enumerate(axes):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylim((0,6.5))
            ax.set_xlabel('absolute change in log2-, z-transf. expression')
            ax.set_ylabel('density')
            ax.annotate('>', (.93, -.06), xycoords='axes fraction')

            m_label = 'measured combinations'
            g_label = 'generated combinations'
            d_label = 'discrepancy'
            e_label = 'explained effect'

            if not i:
                title = ('Distribution of differential effects between measured combinations, \n'
                        'generated combinations and single inhibitors. Combination methods: mean and sum'
                        '\n\ncomined single inhibitor effecs using: `mean`')

                m_data = m_eff.abs().values.flatten()
                g_data = g_eff.abs().values.flatten()
                d_data = n_expl.values.flatten()
                e_data = expl.values.flatten()

                ax.text(1.6, 7, ('genes detected either in generated or measured: 15,640\n'
                        'genes not detected in measured: 1616\n'
                        'genes not detected in generated 2347'),
                        bbox=dict(boxstyle='round', facecolor=config.GREY3)
                )

                
            elif i:
                title = 'comined single inhibitor effecs using: `sum`'
                

                g_ef = self._make_generated_equivalent(method='sum')._expr_eff
                g_ef = g_ef.reindex(genes).drop('Twoi', axis=1, level=1).xs('z', 1, 2)['up']
                g_ef.fillna(0, inplace=True)
                
                n_exp = abs(g_ef - m_ef.values)
                exp = (abs(m_ef.values) - n_exp)
                expl = exp.mask(exp < 0, 0)
                
                n_expl = n_exp.mask(n_exp>2, 2)
                expl = expl.mask(expl>2, 2)
                g_eff = g_ef.mask(g_ef>2, 2).mask(g_ef<-2, -2)
                
                
                m_data = m_eff.abs().values.flatten()
                g_data = g_eff.abs().values.flatten()
                d_data = n_expl.values.flatten()
                e_data = expl.values.flatten()

            ax.set_title(title)

            m_mean = m_data.mean().round(3)
            g_mean = g_data.mean().round(3)
            d_mean = d_data.mean().round(3)
            e_mean = e_data.mean().round(3)

            ax.vlines([m_mean, g_mean, d_mean, e_mean], ymin=0, ymax=4, 
                      color=['purple', '#ffe602', '#ad3232', '#3caf36'], alpha=.6, linewidth=[1.4, 1.4, .4, .4] )
            ax.annotate('mean '+str(m_mean), (m_mean, 3.8), xycoords='data')
            ax.annotate(str(g_mean), (g_mean, 3.2), xycoords='data')
            ax.annotate(str(d_mean), (d_mean, 3), xycoords='data')
            ax.annotate(str(e_mean), (e_mean, 3.8), xycoords='data')
            

            m_max = ax.hist(m_data, density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled', range=(0 , 2), label=m_label)
            g_max = ax.hist(g_data, density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled', range=(0, 2), label=g_label)
            ax.hist(d_data, density=True, bins=200, alpha=.9, color='#ad3232', histtype='step', range=(0, 2), label=d_label)
            e_max = ax.hist(e_data, density=True, bins=200, alpha=.9, color='#3caf36', histtype='step', range=(0, 2), label=e_label)
            print(e_max[0][0])
            
            ax.annotate('explained max:{}'.format(round(e_max[0][0], 1)), xy=(.031, 6.1), xycoords='data')
            ax.annotate('generated max:{}'.format(round(g_max[0][0], 1)), xy=(.031, 5.8), xycoords='data')
            ax.annotate('measured max: {}'.format(round(m_max[0][0], 1)), xy=(.031, 5.5), xycoords='data')
            if i:
                ax.legend(loc='best')

        util.save_file(fig, 'combination discrepancy distributions.png', None)
        plt.show()


extracted effects histogram
        m_ef, g_ef, _, expl, n_exp = data
        m_eff = m_ef.mask(m_ef>2, 2).mask(m_ef<-2, -2)
        g_eff = g_ef.mask(g_ef>2, 2).mask(g_ef<-2, -2)
        expl = expl.mask(expl>2, 2)
        n_expl = n_exp.mask(n_exp>2, 2)

        fig, axes = plt.subplots(nrows=2, figsize=(4, 3))
        fig.subplots_adjust(hspace=.5)

        for i, ax in enumerate(axes):
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.set_ylim((0,6.5))
            ax.set_xlabel('absolute change in log2-, z-transf. expression')
            ax.set_ylabel('density')
            ax.annotate('>', (.93, -.06), xycoords='axes fraction')

            m_label = 'measured combinations'
            g_label = 'generated combinations'
            d_label = 'discrepancy'
            e_label = 'explained effect'

            if not i:
                title = ('Comparison between sum-combined inhibitors and\n'
                         'sum-combined ones + ctracted cross driver effects (3 inhibitor combs only)\n'
                         '\n\nsum-combined inhibitors')

                m_data = m_eff.abs().values.flatten()
                g_data = g_eff.abs().values.flatten()
                d_data = n_expl.values.flatten()
                e_data = expl.values.flatten()

                ax.text(1.6, 7, ('genes detected either in generated or measured: 15,640\n'
                        'genes not detected in measured: 1616\n'
                        'genes not detected in generated 2347'),
                        bbox=dict(boxstyle='round', facecolor=config.GREY3)
                )

                
            elif i:
                title = 'sum-combined inhibitors + extracteted driver effects'

                # c_ef = pd.read_csv('interacs_un.tsv', sep='\t', index_col=0).reindex(g_eff.index)
                
                # cross_eff = pd.read_csv('interacs_notlazy.tsv', sep='\t', index_col=0).reindex(g_eff.index)
                # c_efs = []
                # for comp in self._get_compos():
                #     pce = ['-'.join(sorted([st, nd])) for st in comp 
                #            for nd in comp[comp.index(st)+1:]]
                #     c_efs.append(cross_eff[pce].sum(1).rename('-'.join(comp)))
                # c_ef = pd.concat(c_efs, axis=1)
                

                empt = c_ef.index[c_ef.isna().all(1).values][::2]
                c_ef.fillna(0, inplace=True)
              
                m_ef = m_ef.drop(empt)
                g_ef = g_ef.drop(empt)
                c_ef = c_ef.drop(empt)
                n_exp = n_exp.drop(empt)
                print(g_ef)
                print(c_ef)
                g_ef = c_ef + g_ef.values




                
                n_exp = abs(g_ef - m_ef.values)
                exp = (abs(m_ef.values) - n_exp)
                expl = exp.mask(exp < 0, 0)
                
                n_expl = n_exp.mask(n_exp>2, 2)
                expl = expl.mask(expl>2, 2)
                g_eff = g_ef.mask(g_ef>2, 2).mask(g_ef<-2, -2)
                
                
                m_data = m_eff.abs().values.flatten()
                g_data = g_eff.abs().values.flatten()
                d_data = n_expl.values.flatten()
                e_data = expl.values.flatten()

            ax.set_title(title)

            m_mean = m_data.mean().round(3)
            g_mean = g_data.mean().round(3)
            d_mean = d_data.mean().round(3)
            e_mean = e_data.mean().round(3)

            ax.vlines([m_mean, g_mean, d_mean, e_mean], ymin=0, ymax=4, 
                      color=['purple', '#ffe602', '#ad3232', '#3caf36'], alpha=.6, linewidth=[1.4, 1.4, .4, .4] )
            ax.annotate('mean '+str(m_mean), (m_mean, 3.8), xycoords='data')
            ax.annotate(str(g_mean), (g_mean, 3.2), xycoords='data')
            ax.annotate(str(d_mean), (d_mean, 3), xycoords='data')
            ax.annotate(str(e_mean), (e_mean, 3.8), xycoords='data')
            

            m_max = ax.hist(m_data, density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled', range=(0 , 2), label=m_label)
            g_max = ax.hist(g_data, density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled', range=(0, 2), label=g_label)
            ax.hist(d_data, density=True, bins=200, alpha=.9, color='#ad3232', histtype='step', range=(0, 2), label=d_label)
            e_max = ax.hist(e_data, density=True, bins=200, alpha=.9, color='#3caf36', histtype='step', range=(0, 2), label=e_label)
            print(e_max[0][0])
            
            ax.annotate('explained max:{}'.format(round(e_max[0][0], 1)), xy=(.031, 6.1), xycoords='data')
            ax.annotate('generated max:{}'.format(round(g_max[0][0], 1)), xy=(.031, 5.8), xycoords='data')
            ax.annotate('measured max: {}'.format(round(m_max[0][0], 1)), xy=(.031, 5.5), xycoords='data')
            if i:
                ax.legend(loc='best')

        util.save_file(fig, 'discrepancy cross-effect distri.png', None)

        plt.show()






NAR data discrepany

        m_ef, g_ef, s_ef, expl, n_exp = data
        m_eff = m_ef.mask(m_ef>2, 2).abs()
        s_eff = s_ef.mask(s_ef>2, 2).abs()
        g_eff = g_ef.mask(g_ef>2, 2).abs()
        expl = expl.mask(expl>2, 2)
        n_expl = n_exp.mask(n_exp>2, 2)
        
        # what = m_eff.mask(m_eff<7).dropna(axis=0, how='all')
        # print(m_eff.mask(m_eff<7))
        # return
        fig, axes = plt.subplots(nrows=3, figsize=(6,3))
        fig.subplots_adjust(hspace=.2)
        for i, ax in enumerate(axes):
            if not i:
                title = ('discrepancy distribution of NAR markergenes pooled from all genes of sum-combined inhs. (3 inhs. only)\n\n\n'
                        'all genes (15640) detected in either generated (descending from single inhibitors) or measured')
                ax.set_title(title, pad=8)

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlabel('change in log2-, z-transf. expression')
                ax.set_ylabel('density')
                ax.set_ylim((0,6.5))

                m_mean = round(m_eff.mean().mean(), 3)
                g_mean = round(g_eff.mean().mean(), 3)
                d_mean = round(n_expl.mean().mean(), 3)
                ax.vlines([m_mean, g_mean, d_mean], ymin=(0,0,0), ymax=(3, 2.4, 2),
                color=['purple', '#ffe602', '#ad3232'], alpha=.6, linewidth=(1.4, 1.4, .4))
                ax.annotate('mean '+str(m_mean), (m_mean, 3), xycoords='data')
                ax.annotate('mean '+str(g_mean), (g_mean, 2.4), xycoords='data')
                ax.annotate('mean '+str(d_mean), (d_mean, 2), xycoords='data')

                m_max = ax.hist(m_eff.values.flatten(), density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled', range=(0, 2), label='measured combinations')
                g_max = ax.hist(g_eff.values.flatten(), density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled', range=(0, 2), label='generated combinations')
                ax.hist(n_expl.values.flatten(), density=True, bins=200, alpha=1, color='#ad3232', histtype='step', range=(0, 2), label='discrepancy')
                ax.legend(loc='best')


                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlabel('change in log2-, z-transf. expression', labelpad=-.2)
                ax.set_ylabel('desnsity')

                ax.annotate('generated max:{}'.format(round(g_max[0].max(), 1)), xy=(.031, 6.1), xycoords='data')
                ax.annotate('measured max:{}'.format(round(m_max[0].max(), 1)), xy=(.031, 5.8), xycoords='data')

            elif i ==1:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlabel('change in log2-, z-transf. expression')
                ax.set_ylabel('density')
                ax.set_ylim((0,6.5))

                ax.set_title('NAR markergenes only (8500)', pad=-9)

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.set_xlabel('change in log2-, z-transf. expression', labelpad=-.2)
                ax.set_ylabel('desnsity')


                m_eff = m_eff.reindex(trg._diff_mgs.index)
                g_eff = g_eff.reindex(trg._diff_mgs.index)
                n_expl = n_expl.reindex(trg._diff_mgs.index)

                m_mean = round(m_eff.mean().mean(), 3)
                g_mean = round(g_eff.mean().mean(), 3)
                d_mean = round(n_expl.mean().mean(), 3)

                ax.vlines([m_mean, g_mean, d_mean], ymin=(0,0,0), ymax=(3, 2.4, 2),
                color=['purple', '#ffe602', '#ad3232'], alpha=.6, linewidth=(1.4, 1.4, .4))
                ax.annotate('mean '+str(m_mean), (m_mean, 3), xycoords='data')
                ax.annotate('mean '+str(g_mean), (g_mean, 2.4), xycoords='data')
                ax.annotate('mean '+str(d_mean), (d_mean, 2), xycoords='data')

                m_max = ax.hist(m_eff.values.flatten(), density=True, bins=200, alpha=.6, color='purple', histtype='stepfilled', range=(0, 2), label='measured combinations, NAR genes')
                g_max = ax.hist(g_eff.values.flatten(), density=True, bins=200, alpha=.6, color='#ffe602', histtype='stepfilled', range=(0, 2), label='generated combinations, NAR genes')
                ax.hist(n_expl.values.flatten(), density=True, bins=200, alpha=1, color='#ad3232', histtype='step', range=(0, 2), label='discrepancy, NAR genes')

                ax.annotate('generated max:{}'.format(round(g_max[0].max(), 1)), xy=(.031, 6.1), xycoords='data')
                ax.annotate('measured max:{}'.format(round(m_max[0].max(), 1)), xy=(.031, 5.8), xycoords='data')

            elif i == 2:
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)

                valid = ['all', 'blood mesoderm', 'embryonic', 'endoderm', 'germ cells',
                         'mesoderm', 'neural crest', 'neuroectoderm', 'surface ectoderm']
                n_expls = []
                for i, v in enumerate(valid):

                    t = TARGET(v)
                    lin_mgs = t._diff_mgs.index.intersection(genes)
                    n_expl = n_expl.reindex(lin_mgs)
                    m_n_expl = n_expl.mean().mean()
                    n_expls.append(m_n_expl)
                    ax.annotate(str(round(m_n_expl, 3)), (i, m_n_expl+.05), xycoords='data')
                
                ax.bar(np.arange(len(n_expls)), n_expls, linewidth=1, width=.35, color='#b25555', edgecolor='#ad3232', label='mean discrepancy per lineage-markergene set')


                ax.set_xticks(np.arange(len(valid)))
                ax.set_xticklabels(valid)
                ax.set_xticklabels(valid)
                ax.set_ylim((0, 1))
                ax.set_xlabel('discrepancy')
                ax.legend(loc='upper left')

        plt.show()
        fig.savefig('NAR data discrepancy.png')




        correlation

        plt.xlabel('log2-, z-difference between measured control and generated (single inhs.)\nfor each C/T', fontsize=4)
        plt.ylabel('mean discrepancy across measured/generated combinations', fontsize=4)
        plt.title('Effect of control misalignment between measured combs (-2iLif)\nand generated combs (+2iLif) on discrepancy', fontsize=4)
        
        corr = pd.read_csv('corr_acc_ctrldist_noemb.tsv', sep='\t', index_col=0)
        corr_m = corr.drop('ctrl_dist', axis=1).mean(1)
        ctrl = corr.iloc(1)[-1]
        
        linreg = stats.linregress(ctrl, corr_m)
        print(linreg)
        m = linreg[0]
        b = linreg[1]
        
        plt.plot([0, 2.5], [b, b+(2.5*m)], color='orange', linewidth=.6)
        plt.annotate('r^2 (all, minus embryonic) = {}'.format(round(linreg[2]**2, 4)), (2.3, b+(2.5*m)), fontsize=4)
        # for col in corr.columns[:-1]:
        plt.scatter(ctrl, corr_m, s=1.2, color='orange')
        
        
        
        corrr = pd.read_csv('corr_acc_ctrldist_emb.tsv', sep='\t', index_col=0)
        corr_m = corrr.drop('ctrl_dist', axis=1).mean(1)
        ctrl = corrr.iloc(1)[-1]
        linreg = stats.linregress(ctrl, corr_m)
        m = linreg[0]
        b = linreg[1]
        plt.plot([0, 1.5], [b, b+(1.5*m)], color='green', linewidth=.6)
        plt.annotate('r^2 (embryonic) = {}'.format(round(linreg[2]**2, 4)), (1.3, b+(1.5*m)), fontsize=4)
        plt.scatter(ctrl, corr_m, s=1.2, color='green')
        
        alll = corr.append(corrr)
        
        corr_m = alll.drop('ctrl_dist', axis=1).mean(1)
        ctrl = alll.iloc(1)[-1]
        linreg = stats.linregress(ctrl, corr_m)
        m = linreg[0]
        b = linreg[1]
        plt.plot([0, 2.5], [b, b+(2.5*m)], color='k', linewidth=.6)
        plt.annotate('r^2 (all lineages) = {}'.format(round(linreg[2]**2, 4)), (2.3, b+(2.5*m)), fontsize=4)
        
        # print(alll, corr_m, ctrl)
        
        plt.show()


'''

def obj_fun(ce, *disc_args):
# def obj_fun(ce, d_DMV, d_DMU, d_DMB):
    return (ce[0] + ce[1] + ce[2]) - d_DMV