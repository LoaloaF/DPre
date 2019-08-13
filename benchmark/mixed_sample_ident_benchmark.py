import pandas as pd
import numpy as np
import os
import DPre

os.chdir('benchmark')

"""benchmark multi-target identification with mixed samples (data for F and G)

Loop through passed combinations for mouse and human, generate a dataset 
containing the passed elements in specific proportions. Pass this artifical 
dataset to DPre and score similarity. Return the needeed proportion of the 
original target to reach the threshhold siilarity rank auf perc_thr, defaulting 
to the top 5%.
Writes a TSV for each metric, containing the scores of all combinations.
"""
def benchmark_multi(m_expr, h_expr, m_combs, h_combs, metric, return_single=False, 
                    perc_thr=.05):

    # function that takes the merged expression and combination to synthesize
    # returns the syntheszed expression table
    def synthesize_from_comb(comb, expr):
        # reindex to combination and append a proportion-appendix to the column names
        # that represents how much sample-of-interest is in the synth. sample
        c_exp = expr.reindex(comb.values, axis=1).dropna().fillna(0)
        c_exp.columns = [c+'_0.00' if i != 1 else c+'_1.00' for i, c in enumerate(c_exp.columns)]

        # multiply the single elements with the respective proportions
        base = c_exp.iloc(1)[0].values
        trg = c_exp.iloc(1)[1].values
        if c_exp.shape[1] == 2:
            synth_data = np.array([trg *trg_props[i] + 
                                    base *trg_props[-1-i] 
                                    for i in range(bins)]).T
        else:
            smp_two = c_exp.iloc(1)[2].values
            if c_exp.shape[1] == 3:
                synth_data = np.array([trg *trg_props[i] + 
                                        smp_two *comb2_trg2_props[i] + 
                                        base *comb_base23_props[i]
                                        for i in range(bins)]).T
            elif c_exp.shape[1] == 4:
                smp_three = c_exp.iloc(1)[2].values
                synth_data = np.array([trg *trg_props[i] + 
                                        smp_two *comb3_trg23_props[i] + 
                                        smp_three *comb3_trg23_props[i] + 
                                        base *comb_base23_props[i]
                                        for i in range(bins)]).T
        # set new column names
        cols = ['{}_{:.2f}'.format(c_exp.columns[1][:-5], trg_props[i], 2)
                for i in range(bins)]
        syn_data = pd.DataFrame(synth_data, index=c_exp.index, columns=cols)
        syn_data = pd.concat([c_exp, syn_data], axis=1)
        cols = syn_data.columns.tolist()
        return syn_data.reindex([cols[0]] +cols[2:] +[cols[1]], axis=1)


    res = dict()
    for species, expr, combs in (['mouse', m_expr, m_combs], ['human', h_expr, h_combs]):
        if expr is None:
            continue
        t = DPre.preset_targets(species, color_legend_filename=None)
        # rank must be smaler then min_rank to score in 5% (14/ 16 in mouse and human)
        min_rank = int(len(t) *perc_thr)
        res[species] = dict()
        for i, comb in combs.iterrows():
            comb.dropna(inplace=True)
            # get the optimal, best matching targets
            trg_cts = trg_map.loc[comb.values].values
            
            # built samples-instance from synthetic data
            syn_data = synthesize_from_comb(comb, expr)
            smp = DPre.samples(expression=syn_data, log=False)

            # return all relevent similarity values of samples and targets to 
            # draw a similarity heatmap as in Figure 1 F
            if return_single:
                res[species].update({'_'.join(comb): (t.slice_elements(trg_cts, log=False), smp)})
                return t, smp, trg_cts

            # get the similrity with all targets
            sim = t._get_similarity(smp, metric, log=False)[0].xs('mean', 1, 0)
            asc = True if metric == 'euclid' else False
            # sort the targets to which the combination-synthetic sample is most 
            # similar to. Then get the position of the closest matching target.
            ranks = sim.apply(lambda row: row.sort_values(ascending=asc).index.get_loc(trg_cts[1]), axis=1)
            # req_prop is the synthetic sample with the mininum proportion of 
            # sample-of-interest that was still in the top 5%
            if ranks[-1] <min_rank:
                req_prop = ranks[ranks<min_rank].index[0][-4:]
            else:
                # when higher 5% (i.e. >14 or 16, return the actual rank)
                req_prop = ranks.iloc[-1]
            
            # built return dictionary
            comb_req = pd.Series([req_prop, trg_cts[1]], ['required_prop', 'ref_target'])
            res[species].update({str(i): comb.append(comb_req).to_frame().T})
            print('{}: {}/{}'.format(species, i+1, len(combs)))

        order = ['start', 'target', 'ref_target', 'sample2', 'sample3', 'required_prop']
        res[species] = pd.concat(list(res[species].values()), sort=False).reset_index()
        res[species] = res[species].reindex(order, axis=1)
    return res

# the target map assigns the closest matching target to each sample 
# (ideally, this target is rated as the most transcriptionally similar) 
trg_map = pd.read_csv('mappings/target_mapping.tsv', sep='\t', 
index_col=0, squeeze=True).dropna()

# define the proportions in which combinations are formed
bins = 19
# props of the target of interest
trg_props = np.arange(.05, 1, .05)
# props of the matching base (control/ null-state, of ESCs) for combs with 3 and 4 cell types (2 and 3 somatic noise)
comb_base23_props = np.append(np.arange(.9, .2, -.1), np.arange(.3, 0, -.025))
# props of the ct representing somatic noise in a 3-ct combination
comb2_trg2_props = np.append(np.arange(.05, .4, .05), np.arange(.3, 0, -.025))
# props of the 2 cts representing somatic noise in a 4-ct combination
comb3_trg23_props = np.append(np.arange(.025, .2, .025), np.arange(.15, 0, -.0125))
# glue together
props_arr = np.stack([comb_base23_props, comb2_trg2_props, trg_props], axis=1)
props_arr = np.vstack( [[1,0,0], [0,1,0], props_arr, [0, 0, 1]])

if __name__ == '__main__':
    print(props_arr)
    print(' Sum of rows:\n', props_arr.sum(axis=1)) # should be 1 for all
    h_expr = ['expr/heart_expression.tsv', 
            'expr/hsliver_expression.tsv', 
            'expr/microglia_macro_corticon_expression.tsv']
    # merge human expr to one expression to pool from when making synthetic combinations
    h_expr = pd.concat([pd.read_csv(expr, sep='\t', index_col='ensg') 
                        for expr in h_expr], axis=1, sort=False)
    m_expr = ['expr/ivd_expression.tsv', 
            'expr/iPSC_note_expression.tsv']
    # merge mouse expr to one expression to pool from when making synthetic combinations
    m_expr = pd.concat([pd.read_csv(expr, sep='\t', index_col='ensg') 
                        for expr in m_expr], axis=1, sort=False)

    # RAM overflow issues when doing everything at once
    # cttypes = 'all',
    cttypes = '2ct', '3ct', '4ct'
    metrics = 'euclid', 'pearson', 'cosine'
    # metrics = 'cosine',
    for metric in metrics:
        ct_data = []
        for ct in cttypes:
            h_combs = pd.read_csv('mappings/human_multi_target_mapping_{}.tsv'.format(ct), sep='\t').dropna(how='all')
            m_combs = pd.read_csv('mappings/mouse_multi_target_mapping_{}.tsv'.format(ct), sep='\t').dropna(how='all')
            print(ct, metric, sep='\n')
            os.makedirs('benchmark_scores', exist_ok=True)
            res = benchmark_multi(m_expr, h_expr, m_combs=m_combs, h_combs=h_combs, metric=metric)
            ct_data.extend([res['human'], res['mouse']])
            pd.concat([res['human'], res['mouse']]).to_csv('benchmark_scores/{}_mixed_sample_benchmark_scores_{}.tsv'.format(metric, ct), sep='\t')
        pd.concat(ct_data[::2]+ ct_data[1::2]).to_csv('benchmark_scores/{}_mixed_sample_benchmark_scores_merged.tsv'.format(metric), sep='\t')