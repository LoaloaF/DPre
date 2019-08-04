"""Frontend implementing an argparse command line interface"""
import os
import sys
import argparse

import DPre

# initiate and return the targets and samples for the plots
def init_trg_smp(args):
               if args['preset_targets'] is not None:
                    t = DPre.preset_targets(args['preset_targets'])
               else:
                    t = DPre.targets(markergenes = args['targets_markergenes'],
                         expression = args['targets_expression'],
                         name = args['targets_name'],
                         ignore_down_mgs = args['ignore_down_mgs'],
                         override_namematcher = args['targets_override_namematcher'], 
                         species = args['species'], 
                         log='from_cmd')
               if args['targets_slice'] is not None:
                    t = t.slice_elements(args['targets_slice'])

               s = DPre.samples(diff_genes = args['samples_diff_genes'],
                    expression = args['samples_expression'],
                    ctrl = args['control'],
                    name = args['samples_name'],
                    override_namematcher = args['samples_override_namematcher'], 
                    log='from_cmd')
               if args['samples_slice'] is not None:
                    s = s.slice_elements(args['samples_slice'])
               return t, s

# run the target_similarity_heatmap plot with the parsed argumen ts (1/3)
def _do_target_sim(args):
     t, s = init_trg_smp(args)
     t.target_similarity_heatmap(
                    # plot data
                    samples = s, 
                    which = args['which'], 
                    differential = args['absolute'],
                    proportional = args['proportional'], 
                    display_similarity = args['display_similarity'],
                    # data ordering
                    cluster_targets = args['cluster_targets'],
                    cluster_samples = args['cluster_samples'],
                    reorder_to_distance_bar = args['reorder_to_distance_bar'],
                    # general settings
                    pivot = args['pivot'],
                    heatmap_width = args['heatmap_width'],
                    heatmap_height = args['heatmap_height'],
                    heatmap_range = args['heatmap_range'],
                    distance_bar_range = args['distance_bar_range'],
                    specific_target_labels = args['specific_target_labels'],
                    targetlabels_space = args['targetlabels_space'],
                    samplelabels_space = args['samplelabels_space'],
                    targetlabels_size = args['targetlabels_size'],
                    samplelabels_size = args['samplelabels_size'],
                    title = args['title'],
                    # show/hide elements
                    hide_colorbar_legend = args['hide_colorbar_legend'],
                    hide_distance_bar = args['hide_distance_bar'],
                    hide_targetlabels = args['hide_targetlabels'],
                    hide_targets_dendrogram = args['hide_targets_dendrogram'],
                    hide_targets_colorbar = args['hide_targets_colorbar'],
                    hide_samplelabels = args['hide_samplelabels'],
                    show_samples_dendrogram = args['show_samples_dendrogram'],
                    show_samples_colorbar = args['show_samples_colorbar'],
                    filename = args['filename'])

# run the single gene similarity heatmap with the parsed arguments (3/3)
def _do_gene_sim(args):
     t, s = init_trg_smp(args)
     t.gene_similarity_heatmap(
                    # plot data
                    samples = s,  
                    which = args['which'],
                    differential = args['absolute'],
                    proportional = args['proportional'],
                    display_genes = args['display_genes'],
                    gene_number = args['gene_number'],
                    specific_genes = args['specific_genes'],
                    custom_target_genelist = args['custom_target_genelist'],
                    # data ordering
                    cluster_genes = args['cluster_genes'],
                    cluster_samples = args['cluster_samples'],
                    reorder_to_distance_bar = args['reorder_to_distance_bar'],
                    # general settings
                    pivot = args['pivot'],
                    heatmap_width = args['heatmap_width'],
                    heatmap_height = args['heatmap_height'],
                    heatmap_range = args['heatmap_range'],
                    distance_bar_range = args['distance_bar_range'],
                    sum_plot_range = args['sum_plot_range'],
                    genelabels_space = args['genelabels_space'],
                    samplelabels_space = args['samplelabels_space'],
                    genelabels_size = args['genelabels_size'],
                    title = args['title'],
                    # show/ hide elements
                    hide_colorbar_legend = args['hide_colorbar_legend'],
                    hide_distance_bar = args['hide_distance_bar'],
                    hide_sum_plot = args['hide_sum_plot'],
                    hide_genelabels = args['hide_genelabels'],
                    hide_genes_dendrogram = args['hide_genes_dendrogram'],
                    show_genes_colorbar = args['show_genes_colorbar'],
                    hide_samplelabels = args['hide_samplelabels'],
                    show_samples_dendrogram = args['show_samples_dendrogram'],
                    show_samples_colorbar = args['show_samples_colorbar'],
                    filename = args['filename'])

# run the ranked similarity barplot with the parsed arguments (2/3)
def _do_ranked_sim(args):
     t, s = init_trg_smp(args)
     t.ranked_similarity_barplot(
                    # plot data
                    samples = s,
                    which = args['which'],
                    differential = args['absolute'],
                    proportional = args['proportional'],
                    display_similarity = args['display_similarity'],
                    n_targets = args['n_targets'],
                    display_negative = args['display_negative'],
                    # data ordering
                    rank_samples = args['rank_samples'],
                    # general settings
                    pivot = args['pivot'],
                    xlim_range = args['xlim_range'],
                    targetlabels_space = args['targetlabels_space'],
                    targetlabels_size = args['targetlabels_size'],
                    colored_bars = args['colored_bars'],
                    title = args['title'],
                    spines = args['spines'],
                    # show/ hide elements
                    hide_colorbar = args['hide_colorbar'],
                    hide_targetlabels = args['hide_targetlabels'],
                    filename = args['filename'])

# create the base parser
d = ('DPre - visualizing transcriptional similarity between samples and targets')
u = ('\n\n\n\n\n\nPass the input --> Choose the plot --> specify the plot:\n'
    '$ ./dpre.py <targets> <samples> <plot> <plot args> (optional)\nExample:\n'
     'python dpre.py -pt "h endoderm" -se "examples/example1_hsliver/hsliver_e'
     'xpression.tsv" -c "Day00" target_sim -w "euclid"')
parser = argparse.ArgumentParser(description=d, usage=u, allow_abbrev=False,
                                 add_help=False)

# create a group from the parser for the target data input
trg_grp = parser.add_argument_group('targets', description='input comparison '
                                    'data, pass `preset_targets` or `targets_'
                                    'markergenes` & `targets_expression`')
trg_grp.add_argument('--preset_targets', '-pt', type=str, 
                     help='load a default targets profile (all other options '
                     'are ignored)')
trg_grp.add_argument('--targets_markergenes', '-tm', action='append', type=str, 
                     help='directory/ies with markergene files (up or up & '
                     'down)')
trg_grp.add_argument('--targets_expression', '-te', type=str, 
                     help='filename of targets tsv expression table')
trg_grp.add_argument('--ignore_down_mgs', '-i', action='store_true', 
                     help='even if found in input, don`t use down markergenes')
trg_grp.add_argument('--targets_name', '-tn', type=str, 
                     help='targets name used in headers and logs')
trg_grp.add_argument('--targets_override_namematcher', '-to', 
                     action='store_true', help='when both markergenes and '
                     'expression, ignore name mismatches')
trg_grp.add_argument('--species', '-s', type=str, default='mouse',
                     help='species of the target, default `mouse`')
trg_grp.add_argument('--targets_slice', '-ts', action='append', type=str, 
                     help='convenience slicer, pass the element names to keep')

# create a group from the parser for the samples data input
smp_grp = parser.add_argument_group('samples', description='input data to '
                                    'explore similarity for, pass `targets_diff'
                                    '_genes` and/or `samples_expression`')
smp_grp.add_argument('--samples_diff_genes', '-sd', action='append', type=str, 
                     help='directory/ies with differential gene files (up or '
                     'up & down)')
smp_grp.add_argument('--samples_expression', '-se', type=str, 
                     help='filename of samples tsv expression table')
smp_grp.add_argument('--control', '-c', type=str, 
                     help='samples control name in data')
smp_grp.add_argument('--samples_name', '-sn', type=str, 
                     help='samples name used in headers and logs')
smp_grp.add_argument('--samples_override_namematcher', '-so', 
                     action='store_true', help='when both differential and '
                    'expression, ignore name mismatches')
smp_grp.add_argument('--samples_slice', '-ss', action='append', type=str, 
                     help='convenience slicer, pass the element names to keep')

# create a subparser of parser that specifies the plot to run
# each plot is a parser of this subparser and implements its do_plot function
d = 'one of 3 plots for similarity visualization, for details run <plot> -h'
subparsers = parser.add_subparsers(title='Plots', description=d)
subparsers.required = True

# add help argument here to show at bottom instead of top
optional = parser.add_argument_group('optional arguments')
optional.add_argument('-h', '--help', action='help', 
                      help='show this help message and exit')



# create the parser of the target_similarity_hm; parse respective args
d = ('Plot the elementwise similarity of the samples and targets in a heatmap\n'
     'target_similarity_heatmap specific args:')
trg_sim_parser = subparsers.add_parser('target_sim', description=d, usage=u)
trg_sim_parser.set_defaults(func=_do_target_sim)
trg_sim_parser.add_argument('--filename', '-f', type=str, default='target_'
                            'similarity_hm.png', help='filename for saving.' 
                            'If None plot is not saved, default: '
                            'target_similarity_hm.png')
# add arguments for each argument group
# plot data
d = 'main parameters to control the presented similarity'
dat_grp = trg_sim_parser.add_argument_group('Data options', description=d)
dat_grp.add_argument('--which', '-w', type=str, choices=('euclid', 'intersect'),
                     help='select the similarity metric')
dat_grp.add_argument('--absolute', '-a', action='store_false', help='plot the '
                     'absolute similarity, `euclid` only. Will appear as '
                     'differential=False in logs.')
dat_grp.add_argument('--proportional', '-p', action='store_true',
                     help='plot the proportional changes in similarity')
dat_grp.add_argument('--display_similarity', '-ds', default='mgs mean',
                     choices=['mgs mean', 'mgs up', 'mgs down'], help='Specify '
                     'up- or down markerene similarity, default mean')
# data ordering
d = 'parameters to control ordering, i.e. clustering'
datord_grp = trg_sim_parser.add_argument_group('Data order options', 
                                               description=d)
datord_grp.add_argument('--cluster_targets', '-ct', action='store_true',
                        help='cluster targets (x axis)')
datord_grp.add_argument('--cluster_samples', '-cs', action='store_true',
                        help='cluster samples (y axis)')
datord_grp.add_argument('--reorder_to_distance_bar', '-re', action='store_true',
                        help=('reorder the targets according to the distance '
                              'bar for specific target identification'))
# general heatmap settings
d = 'parameters to control general visual options'
genhm_grp = trg_sim_parser.add_argument_group('General heatmap options', 
                                              description=d)
genhm_grp.add_argument('--pivot', '-pi', action='store_true',
                       help='flip the plot 90 degrees')
genhm_grp.add_argument('--heatmap_width', '-hw', type=float,
                       help='heatmap width multiplier, deafult 1')
genhm_grp.add_argument('--heatmap_height', '-hh', type=float,
                       help='heatmap height multiplier, deafult 1')
genhm_grp.add_argument('--heatmap_range', '-hr', action='append', type=float,
                       help='range of heatmap values, (lower, upper)')
genhm_grp.add_argument('--distance_bar_range', '-dr', action='append', type=float,
                       help='range of distance bar values, (lower, upper)')
genhm_grp.add_argument('--specific_target_labels', '-st', action='append', 
                       type=str,  help='plot specific set of target labels')
genhm_grp.add_argument('--targetlabels_space', '-ta', type=float, 
                       help='space reserved for targetlabels in inches')
genhm_grp.add_argument('--samplelabels_space', '-sa', type=float, 
                       help='space reserved for sample labels in inches')
genhm_grp.add_argument('--targetlabels_size', '-tas', type=float, 
                       help='multiplier for targetlabels fontsize, default = 1')
genhm_grp.add_argument('--samplelabels_size', '-sas', type=float, 
                       help='multiplier for samplelabels fontsize, default = 1')
genhm_grp.add_argument('--title', '-t', default=True,
                       help='a custom title or hide title if `f`, `F`, ..')
# show/ hide specific plot elements
d = 'show/ hide subparts of the plot'
elem_grp = trg_sim_parser.add_argument_group('Plot elements', description=d)
elem_grp.add_argument('--hide_colorbar_legend', '-hco', action='store_true', 
                      help='do not plot the heatmap legend')
elem_grp.add_argument('--hide_distance_bar', '-hd', action='store_true',
                      help='do not plot the distance bar')
elem_grp.add_argument('--hide_targetlabels', '-hta', action='store_true', 
                      help='do not show the targetlabels')
elem_grp.add_argument('--hide_targets_dendrogram', '-htd', action='store_true',
                      help='do not plot the targets dendrogram')
elem_grp.add_argument('--hide_targets_colorbar', '-htc', action='store_true',
                      help='do not plot the targets colorbar')
elem_grp.add_argument('--hide_samplelabels', '-hsa', action='store_true', 
                      help='do not show the samplelabels')
elem_grp.add_argument('--show_samples_dendrogram', '-ssd', action='store_true',
                      help='do not plot the samples dendrogram')
elem_grp.add_argument('--show_samples_colorbar', '-ssc', action='store_true',
                      help='do not plot the samples colorbar')


# create the parser of the gene_similarity_hm; parse respective args
d = 'Plot the single-gene similarity of the samples with each target in a heatmap'
gene_sim_parser = subparsers.add_parser('gene_sim', description=d, usage=u)
gene_sim_parser.set_defaults(func=_do_gene_sim)
gene_sim_parser.add_argument('--filename', '-f', type=str, default='gene_'
                             'similarity_hm.pdf', help='filename for saving.' 
                             'If None plot is not saved, default: '
                             'gene_similarity_hm.pdf')
# add arguments for each argument group
# plot data
d = 'main parameters to control the presented similarity'
dat_grp = gene_sim_parser.add_argument_group('Data options', description=d)
dat_grp.add_argument('--which', '-w', type=str, choices=('euclid', 'intersect'),
                     help='select the similarity metric')
dat_grp.add_argument('--absolute', '-a', action='store_false', help='plot the '
                     'absolute similarity, `euclid` only. Will appear as '
                     'differential=False in logs.')
dat_grp.add_argument('--proportional', '-p', action='store_true',
                     help='plot the proportional changes in similarity')
dat_grp.add_argument('--display_genes', '-di', type=str, default='variant', 
                     choices=['variant', 'increasing', 'decreasing', 
                     'similar', 'distant'], help='Option 1: specify the group '
                     'of genes to extract, default `variant`')
dat_grp.add_argument('--gene_number', '-gn', type=int, default=45, 
                     help='specify number of genes to extract (for to opt. 1)')
dat_grp.add_argument('--specific_genes', '-sp', action='append', 
                     help='Option 2: specify a list of target markergenes' )
dat_grp.add_argument('--custom_target_genelist', '-cu',  action='append', 
                     help='Option 3: specify any list of genes')
# data ordering
d = 'parameters to control ordering, i.e. clustering'
datord_grp = gene_sim_parser.add_argument_group('Data order options', 
                                                description=d)
datord_grp.add_argument('--cluster_genes', '-cg', action='store_true',
                        help='cluster genes (x axis)')
datord_grp.add_argument('--cluster_samples', '-cs', action='store_true',
                        help='cluster samples (y axis)')
datord_grp.add_argument('--reorder_to_distance_bar', '-re', action='store_true',
                        help=('reorder the genes according to the distance '
                        'bar, i.e the base gene similarity'))
# general heatmap settings
d = 'parameters to control general visual options'
genhm_grp = gene_sim_parser.add_argument_group('General heatmap options', 
                                               description=d)
genhm_grp.add_argument('--pivot', '-pi', action='store_true',
                       help='flip the plot 90 degrees')
genhm_grp.add_argument('--heatmap_width', '-hw', type=float,
                       help='heatmap width multiplier, default 1')
genhm_grp.add_argument('--heatmap_height', '-hh', type=float,
                       help='heatmap height multiplier, default 1')
genhm_grp.add_argument('--heatmap_range', '-hr', action='append', type=float, 
                       help='range of heatmap values, (lower, upper)')
genhm_grp.add_argument('--distance_bar_range', '-dr', action='append', type=float, 
                       help='range of distance bar values, (lower, upper)')
genhm_grp.add_argument('--sum_plot_range', '-sr', action='append', type=float,
                       help='range of sum plot values, (lower, upper)')
genhm_grp.add_argument('--genelabels_space', '-ge', type=float, 
                       help='space reserved for genelabels in inches')
genhm_grp.add_argument('--samplelabels_space', '-sa', type=float, 
                       help='space reserved for samplelabels on left in inches')
genhm_grp.add_argument('--genelabels_size', '-ges', type=float, 
                       help='multiplier for genelabels fontsize, default = 1')
genhm_grp.add_argument('--targetlabels_size', '-tas', type=float, 
                       help='multiplier for targetlabels fontsize, default = 1')
genhm_grp.add_argument('--title', '-t', default=True, help='a custom title or '
                       'hide title with {f, false, F, False}')
# show/ hide specific plot elements
d = 'show/ hide subparts of the plot'
elem_grp = gene_sim_parser.add_argument_group('Plot elements', description=d)
elem_grp.add_argument('--hide_colorbar_legend', '-hco', action='store_true', 
                      help='do not plot the legend for the heatmap')
elem_grp.add_argument('--hide_distance_bar', '-hd', action='store_true',
                      help='do not plot the distance bar')
elem_grp.add_argument('--hide_sum_plot', '-hs', action='store_true',
                      help='do not plot the sum plot on the right')
elem_grp.add_argument('--hide_genelabels', '-hge', action='store_true', 
                      help='do not show the genelabels')
elem_grp.add_argument('--hide_genes_dendrogram', '-hgd', action='store_true',
                       help='do not plot the genes dendrogram')
elem_grp.add_argument('--show_genes_colorbar', '-sgc', type=bool, default=False,
                      help='plot a genes colorbar when specific_genes passed')
elem_grp.add_argument('--hide_samplelabels', '-hsa', action='store_true', 
                      help='do not show the samplelabels')
elem_grp.add_argument('--show_samples_dendrogram', '-ssd', action='store_true',
                      help='do plot the samples dendrogram')
elem_grp.add_argument('--show_samples_colorbar', '-ssc', action='store_true',
                      help='do plot the samples colorbar')



# create the parser of the ranked_similarity_barplot; parse respective args
d = 'Plot the ranked similarity of the samples with targets in a bar plot'
rank_sim_parser = subparsers.add_parser('ranked_sim', description=d, usage=u)
rank_sim_parser.set_defaults(func=_do_ranked_sim)
rank_sim_parser.add_argument('--filename', '-f', type=str, default='ranked_'
                             'similarity_bp.pdf', help='filename for saving.' 
                             'If None plot is not saved, default: '
                             'ranked_similarity_hm.png')
# add arguments for each argument group
# plot data
d ='main parameters to control the presented similarity'
dat_grp = rank_sim_parser.add_argument_group('Data options', description = d)
dat_grp.add_argument('--which', '-w', type=str, choices=('euclid', 'intersect'),
                     help='select the similarity metric')
dat_grp.add_argument('--absolute', '-a', action='store_false', help='plot the '
                     'absolute similarity, `euclid` only. Will appear as '
                     'differential=False in logs.')
dat_grp.add_argument('--proportional', '-p', action='store_true',
                     help='plot the proportional changes in similarity')
dat_grp.add_argument('--display_similarity', '-ds', default='mgs mean',
                     choices=['mgs mean', 'mgs up', 'mgs down'], help='Specify '
                     'up- or down markerene similarity, default mean')
dat_grp.add_argument('--n_targets', '-nt', type=int, default=16,
                     help='specify number of targets to show')
dat_grp.add_argument('--display_negative', '-din', action='store_true',
                     help='besides postive values, display the negative ones')
# data ordering
d = 'parameters to control ordering'
datord_grp = rank_sim_parser.add_argument_group('Data order options', 
                                                description=d)
datord_grp.add_argument('--rank_samples', '-ra', action='store_true', help='re'
                        'order the plots according to the peak sample values')
# general barplot settings
d = 'parameters to control general visual options'
genhm_grp = rank_sim_parser.add_argument_group('General barplot options', 
                                               description=d)
genhm_grp.add_argument('--pivot', '-pi', action='store_true',
                       help='flip the plot 90 degrees')
genhm_grp.add_argument('--xlim_range', '-x', action='append', type=float, 
                       help='range of barplot values, (lower, upper)')
genhm_grp.add_argument('--targetlabels_space', '-ta', type=float, 
                       help='space reserved for targetlabels in inches')
genhm_grp.add_argument('--targetlabels_size', '-tas', type=float, 
                       help='multiplier for targetlabels fontsize, default = 1')
genhm_grp.add_argument('--colored_bars', '-co', action='store_true',
                       help='color the bars according to the value')
genhm_grp.add_argument('--spines', '-spi', action='store_true', 
                      help='additionally show spines on top and right')
genhm_grp.add_argument('--title', '-t', default=True, help='a custom title or '
                       'hide title with {f, false, F, False}')
# show/ hide specific plot elements
d = 'show/ hide subparts of the plot'
elem_grp = rank_sim_parser.add_argument_group('Plot elements', description=d)
elem_grp.add_argument('--hide_targetlabels', '-hta', action='store_true', 
                      help='do not show the target labels')
elem_grp.add_argument('--hide_colorbar', '-hc', action='store_true', 
                      help='do not plot the targets colorbar')

args = vars(parser.parse_args())
do_plot = args.pop('func')
do_plot(args)