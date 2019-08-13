import pandas as pd
import os
import DPre

os.chdir('benchmark')

"""benchmark target identification (Figure 1E) 
    Writes a TSV file that contains the similarity at which the closest 
    matching target for each sample was ranked. Optimal is 1. Runs for all 4 
    metrics. 
"""

def benchmark_single(metrics):
    scores = dict(euclid=[], pearson=[], cosine=[], intersect=[])
    for name_expr, t in species:
        for smp_name, expr, ctrl in name_expr:
            s = DPre.samples(expression=expr, name=smp_name, ctrl=ctrl, log=False)
            if metrics == ['intersect']:
                DPre.add_diff_genes_from_z(s)

            for metric in metrics:
                diff = False if metric != 'intersect' else True
                # get the similarity of all samples with all targets
                sim = t._get_similarity(s, metric, differential=diff, log=False)[0].xs('mean', 1, 0)
                # slice to all relevent targets (the closest matching ones)
                sim = sim.reindex(sim.index.intersection(trg_map.index))
                # iterate the samples
                for name, row in sim.iterrows():
                    asc = True if metric == 'euclid' else False 
                    # get rank position of best matching target
                    # print(row.sort_values(ascending=asc))
                    rank = row.sort_values(ascending=asc).index.get_loc(trg_map[name]) +1
                    print(rank)
                    scores[metric].append([name, trg_map[name], rank])
                    print(metric, ': ', name, ' done.')
            print()
    # assamble a results dataframe
    scores_df = []
    for key in scores:
        res = pd.DataFrame({'sample': [score[1] for score in scores[key]], 
                            'rank': [score[2] for score in scores[key]], 
                            'metric': [key]*len(scores[key])}, 
                            index=[score[0] for score in scores[key]])
        scores_df.append(res)
    return pd.concat(scores_df, axis=0).set_index('metric', append=True).swaplevel()

# the target map assigns the closest matching target to each sample 
# (ideally, this target is rated as the most transcriptionally similar) 
trg_map = pd.read_csv('mappings/target_mapping.tsv', sep='\t', 
                      index_col=0, squeeze=True).dropna()

# non-differenital metrics                      
h_expr = [('Liver', 'expr/hsliver_expression.tsv', None),
          ('Heart' , 'expr/heart_expression.tsv', None), 
          ('MacroMicroNeuro', 'expr/microglia_macro_corticon_expression.tsv', None)]
m_expr = [('IVD', 'expr/ivd_expression.tsv', None),
          ('MEF', 'expr/iPSC_note_expression.tsv', None)]
species = ((h_expr, DPre.preset_targets('human', color_legend_filename=None)), 
           (m_expr, DPre.preset_targets('mouse', color_legend_filename=None)))
# run benchmark that scores the rank of the best matching target reference (rank 1 is optimal)
nondiff = benchmark_single(['euclid', 'cosine', 'pearson'])

# differential metrics (only intersect)
m_expr = [('IVD', 'expr/ivd_expression_mneu.tsv', 'ivd ESCs (mneu)'),
          ('IVD2', 'expr/ivd_expression_cardio.tsv', 'ivd ESCs (cardio)'),
          ('MEF', 'expr/iPSC_note_expression.tsv', 'MEF')]
h_expr = [('Liver', 'expr/hsliver_expression.tsv', 'Day00'),
          ('Heart' , 'expr/heart_expression.tsv', 'mean_hESC'), 
          # because this dataset doesn't have a control, the control from the Heart data was added 
          ('MacroMicroNeuro', 'expr/microglia_macro_corticon_expression_art_ctrl.tsv', 'mean_hESC')]
species = ((h_expr, DPre.preset_targets('human', color_legend_filename=None)), 
           (m_expr, DPre.preset_targets('mouse', color_legend_filename=None)))
# run benchmark that scores the rank of the best matching target reference (rank 1 is optimal)
diff = benchmark_single(['intersect'])

nondiff.append(diff).to_csv('benchmark_scores/ct_ident_benchmark_scores.tsv', sep='\t')