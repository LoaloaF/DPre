from DPre import samples
from DPre import preset_targets
import pandas as pd
from DPre import targets    

# iniatite multiple sample data instances because of different controls
ivd_expr = pd.read_csv('ivd_expression.tsv', sep='\t', index_col='ensg')
# cardiomyocytes 
cm_sample = samples(expression = ivd_expr.loc[:, ['ivd ESCs (cardio)', 'ivd cardiomyocytes']],    # syntax: loc[all rows, specifc columns] 
                    ctrl = 'ivd ESCs (cardio)',
                    name = 'in vitro differentiated cardiomyocytes')
# motor neurons
mn_sample = samples(expression = ivd_expr.loc[:, ['ivd ESCs (mneu)', 'ivd motor neurons']], 
                    ctrl = 'ivd ESCs (mneu)',
                    name = 'in vitro differentiated motor neurons')
# retinoic acid treated
ra_sample = samples(expression = ivd_expr.loc[:, ['ivd ESCs (ra)', 'ivd retinoic acid']],
                    ctrl = 'ivd ESCs (ra)',
                    name = 'in vitro differentiated ESCs +retinoic acid')

# initatie the target to compare against and check target marker gene detection proportion 
t = preset_targets('mouse')
hist = t.plot_detec_mgs_prop(samples = cm_sample,
                             filename = 'trg_mgs_detec_prop.png',
                            )

# make the `in vitro` in the labels italic
it_in_vit = '$\mathit{in}$ $\mathit{vitro}$'     # TEX expression for italic 'in vitro' string
cm_sample.names = ({'ivd cardiomyocytes': it_in_vit+' differentiated cardiomyocytes'})    # passing a mapping of old name -> new name
mn_sample.names = ({'ivd motor neurons': it_in_vit+' differentiated motor neurons'})
ra_sample.names = ({'ivd retinoic acid': 'ESCs +retinoic acid'})

# draw the ranked similarity plot iteratively
for sam in (cm_sample, mn_sample, ra_sample):
    t.ranked_similarity_barplot(samples = sam,
                                metric='euclid', 
                                n_targets = 10, 
                                xlim_range = [-1.7, 1.7],    # ensure consistent ranges across plots
                                display_negative = True,
                                pivot = True,
                                targetlabels_size =.95,      # downscale the labels for space saving
                                targetlabels_space= 1,
                                BP_BARWIDTH_SIZE = .14,
                                BP_BOTTOM = .2,
                                BP_TOP = 1,
                                title = sam.names[1],       # set the element name as the title
                                filename = 'ranked_sim_'+ sam.name +'.svg',    # when .png, DPre autiamically uniquely extends the filename
                                )
# commands
# python ../../dpre.py -pt "mouse" -se "ivd_expre> python ../../dpre.py -pt "mouse" -se "ivd_expression.tsv" -c "ivd ESCs (cardio)" -ss "ivd cardiomyocytes"  -ss "ivd ESCs (cardio)"  -sn "in vitro differentiated cardiomyocytes" ranked_sim -nt 10 -x -1.7 -x 1.7 -din -pi -tas .95 -ta 1 -t "in vitro differentiated cardiomyocytes" -f "ranked_sim_in vitro differentiated cardiomyocytes.svg"ssion.tsv" -c "ivd ESCs (cardio)" -ss "ivd cardiomyocytes"  -ss "ivd ESCs (cardio)"  -sn "in vitro differentiated cardiomyocytes" ranked_sim -nt 10 -x -1.7 -x 1.7 -din -pi -tas .95 -ta .6 -f "ranked_sim.png"
# python ../../dpre.py -pt "mouse" -se "ivd_expression.tsv" -c "ivd ESCs (mneu)" -ss "ivd motor neurons"  -ss "ivd ESCs (mneu)" -sn "in vitro differentiated motor neurons" ranked_sim -nt 10 -x -1.7 -x 1.7 -din -pi -tas .95 -ta 1 -t "in vitro differentiated motor neurons" -f "ranked_sim_vitro differentiated motor neurons.svg"
# python ../../dpre.py -pt "mouse" -se "ivd_expression.tsv" -c "ivd ESCs (ra)" -ss "ivd retinoic acid" -ss "ivd ESCs (ra)" -sn "in vitro differentiated mESCs +retinoic acid" ranked_sim  -nt 10 -x -1.7 -x 1.7 -din -pi -tas .95 -ta 1 -t "in vitro differentiated mESCs +retinoic acid" -f "ranked_sim_in vitro differentiated mESCs +retinoic acid.svg"

# make a samples instance with the 3 samples above to produce one heatmap
all_expr = ivd_expr.loc[:, ['ivd cardiomyocytes', 'ivd motor neurons', 'ivd retinoic acid']]      # exclude controls
all_expr.columns = [it_in_vit+' differentiated cardiomyocytes',     # rename expression input for directly getting the correct names
                    it_in_vit+' differentiated motor neurons', 
                    'ESCs +retinoic acid']                                   
all_sample = samples(expression = all_expr, 
                     name = it_in_vit + ' differentiated ESCs')
hm = t.target_similarity_heatmap(samples = all_sample,           # because no metric is passed, the default metric for target similarity is used: consine similarty
                                 differential = False,             # on command line this is passed by absolute = True, i.e. -a
                                 heatmap_height = 1.6,        # make heatmap higher
                                 heatmap_width = .13, 
                                 hide_targetlabels = True, 
                                 targetlabels_space = .45,
                                 samplelabels_space = 1.7,
                                 HM_TOP = 1,
                                 filename = 'abs_sim_hm.png',
                                 )
# command
# python ../../dpre.py -pt "mouse" -se "ivd_expression.tsv" -ss "ivd cardiomyocytes"  -ss "ivd motor neurons"  -ss "ivd retinoic acid" -sn "in vitro differentiated ESCs" target_sim -a -hh 1.6 -hw .13 -hta -ta .45 -sa 1.7 -f "abs_sim_hm.png"