'''

Clean up the human FANTOM5 data.

This version is the non-cancer cell containing one.

You will need glbase3 to run this:

https://bitbucket.org/oaxiom/glbase3

'''

import gzip
from glbase3 import *

# This is the CAGE peak annotations:
oh = gzip.open('hg38_fair+new_CAGE_peaks_phase1and2_ann.txt.gz', 'rt') # Downloaded from FANTOM5, see the accompanying bash script

gene_annots = []
for idx, line in enumerate(oh):
    if '#' in line[0]:
        continue

    line = line.strip().split('\t')
    #print(line)

    if '@' in line[1] and 'ENST0' in line[3]: # In this data we only care about the gene annotated CAGE:
        ent = {'fantom_name': line[0].split(';')[1],
            'promoter_usage': line[1],
            'enst': line[3].split('_')[2].split('.')[0]}
        gene_annots.append(ent)

    if (idx+1) % 100000 == 0:
        print('Processed: {:,} annotations'.format(idx+1))
        #break
oh.close()
gene_annots = genelist(loadable_list=gene_annots)
valid_ids = set(gene_annots['fantom_name'])

print('Found {:,} valid annotations'.format(len(gene_annots)))

oh = gzip.open("hg38_fair+new_CAGE_peaks_phase1and2_tpm.osc.txt.gz", "rt")# Downloaded from FANTOM5, see the accompanying bash script

res = []
cond_names = []
reps = {}

for entry, line in enumerate(oh):
    if "00Annotation" in line: # get the condition names
        tt = line.strip().split("\t")
        cond_names = ["".join(i.split(".")[0:2]) for i in tt[7:]]
        cond_names = [i.replace("%20", " ").replace("tpm", "").replace("%2b", "+").replace("%2c", "").replace("%3a", ":") for i in cond_names]
        cond_names = [i.replace("%28", "(").replace("%29", ")").replace("%27", "").replace('%2f', '-') for i in cond_names]

        # And also remove time course data, (response to), controls (Univeral RNA) and others:
        cancer_names = set([
            "lymphoma", "carcinoma", "cancer", "infection", # cancers;
            "lipopolysaccharide", "Hep-2", "neuroblastoma", "tumor", "sarcoma", "SABiosciences",
            "treated",  "myotonica", "atrophy", "walker warburg", "Embryoid_body",
            "cell line", 'differentiation', 'periodontitis', 'differentiated', 'post-infarction', 'asthmatic',
            "nuclear fraction", "Universal RNA", 'Whole blood',"Clontech", # controls;
            "induction", 'response to', 'induced with', # time course treatments;
            ])
        to_keep = [cond_names.index(i) for i in cond_names if True not in [t in i for t in cancer_names]]
        cond_names = [cond_names[idx] for idx in to_keep]
        continue

    if (entry+1) % 10000 == 0:
        print('Processed: {:,} genes'.format(entry+1))
        #break

    if '#' in line:
        continue

    if 'hg_' not in line:
        continue

    tt = line.strip().split("\t")
    fantom_name = tt[0].split(";")[1] # the first part is the old name, the new part is the old name;
    if fantom_name not in valid_ids:
        continue

    res.append({"fantom_name": fantom_name, "conditions": [float(tt[7:][idx]) for idx in to_keep]})

# Keep only the genes with an annotated GENE ENST:
gl = expression(loadable_list=res, cond_names=cond_names)
gl = gene_annots.map(genelist=gl, key='fantom_name')

cuts = ['donor', 'biol_rep', 'tech_rep', '(', ')', 'rep', 'donation', 'pool', 'biol_', 'and', '-']
for k in cond_names:
    head = []
    for t in k.split(' '):
        if True in [c in t for c in cuts]:
            continue
        head.append(t)
    head = ' '.join(head).capitalize()

    if head not in reps:
        reps[head] = []
    reps[head].append(k)
#print(reps)

# mean replicates
gl = gl.mean_replicates(*reps.values(),
    threshold=0.6,
    output_pears="Pearson_correlations.tsv",
    pearson_hist="Pearson_histogram.png")
gl.strip_errs()

# Do some sample name tidying:
newconds = []
for c in gl.getConditionNames():
    newc = []
    for t in c.split(' '):
        if True in [c in t for c in cuts]:
            continue
        newc.append(t)
    newc = ' '.join(newc)

    if newc in newconds:
        print('ERROR: clashing condition names: "%s" vs "%s"' % (c, newc))

    newconds.append(newc.capitalize())

gl.setConditionNames(newconds)

gl.sort_sum_expression() # resort the data so picking will get the highest expressed promoter
gl.reverse()
gl.save("hg38.fantom5.glb")
gl.saveTSV("hg38.fantom5.tsv")

print("Number of tss: {:,}".format(len(gl)))

# Some diagnostics:
gl.tree(filename="tree.png", color_threshold=0.0, label_size=4, size=(5,18))
gl.correlation_heatmap(filename="corr_heatmap.png", bracket=(0.0,1.0), size=(17,16), heat_wid=0.64, heat_hei=0.927, row_font_size=5)
gl.network.conditions(filename="ctype_net.png", low_threshold=0.60, hi_threshold=0.80,
    max_links=10, size=(16,16))
