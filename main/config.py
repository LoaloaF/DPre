"""Module storing package wide variables, colors and plot sizes. May be edited
by the user"""
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import set_link_color_palette

# change DPI from 100 to 300, adjust everything accordingly
mpl.rcParams['figure.dpi'] = 300.0
mpl.rcParams['savefig.dpi'] = 300.0
mpl.rcParams['font.size'] = FONTS = 5
mpl.rcParams['lines.linewidth'] = 0.3
mpl.rcParams['hatch.linewidth'] = 0.3333
mpl.rcParams['patch.linewidth'] = 0.3333
mpl.rcParams['axes.linewidth'] = 0.26664
mpl.rcParams['xtick.major.size'] = 1.16655
mpl.rcParams['xtick.major.width'] = 0.26664
mpl.rcParams['ytick.major.size'] = 1.16655
mpl.rcParams['ytick.major.width'] = 0.26664
mpl.rcParams['grid.linewidth'] = 0.26664
mpl.rcParams['xtick.major.pad'] = 1.16655
mpl.rcParams['ytick.major.pad'] = 1.16655
mpl.rcParams['figure.max_open_warning'] = 30

# the default save format if the passed filename has no valid ending
SAVE_FORMAT = 'pdf'

# To give relevant insight on transcriptional similarity, a minimum 
# proportion of the targets markergenes should be detected in the samples.
# If the markergene detection proportion is below the specified value, this
# target is dropped from the analysis. By default 15%.
DROP_TARGET_DETEC_THR = .15

# p value when dealing with deseq2 gene list input
DESEQ2_P = .05

# predefined colors to use for labeling 
colors = ['#e6194B', #    0 = red
          '#3cb44b', #    1 = green
          '#ffe119', #    2 = yellow
          '#4363d8', #    3 = blue
          '#f58231', #    4 = orange
          '#911eb4', #    5 = purple
          '#42d4f4', #    6 = cyan
          '#f032e6', #    7 = magenta
          '#bfef45', #    8 = lime
          '#fabebe', #    9 = pink
          '#469990', #    10 = teal
          '#e6beff', #    11 = lavender
          '#9A6324', #    12 = brown
          '#fffac8', #    13 = beige
          '#8a0b25', #    14 = deep red
          '#aaffc3', #    15 = mint
          '#808000', #    16 = olive
          '#ffd8b1', #    17 = apricot
          '#0a3b70', #    18 = deep blue
          '#a9a9a9', #    19 = grey
          '#ffffff', #    20 = white
          '#000000'  #    21 = black
]

# predefined colors for preset_targets
preset_targets_colors = {
    'm embryonic': '#fbcc8c', 
    'm germ cells': '#f0f0a6',
    'm neural crest': '#c8b1ce', 
    'm surface ectoderm': '#9edce4',
    'm neuroectoderm': '#5888c2', 
    'm mesoderm': '#64bb79', 
    'm endoderm': '#843487', 
    'm blood mesoderm': '#fe7e81', 
    'h embryonic': '#fbcc8c', 
    'h surface ectoderm': '#9edce4',
    'h neuroectoderm': '#5888c2', 
    'h mesoderm': '#64bb79', 
    'h endoderm': '#843487', 
    'h blood mesoderm': '#fe7e81', 
}

# predefined legend parameters for annotating preset_targets with various color-
# labeled groups. For presets below, colors are defined in colors.tsv.
preset_col_legend = {
    'mouse': (['Embryonic', 'Germ cells', 'Neural crest', 'Surface ectoderm', 
               'Neuroectoderm', 'Mesoderm', 'Endoderm', 'Blood mesoderm'], 
              ['#fbcc8c', '#f0f0a6', '#c8b1ce', '#9edce4', '#5888c2',  
               '#64bb79',  '#843487',  '#fe7e81']),
    'human': (['Embryonic', 'Surface ectoderm', 'Neuroectoderm', 'Mesoderm',  
               'Endoderm', 'Blood mesoderm'],
              ['#fbcc8c', '#9edce4', '#5888c2', '#64bb79', '#843487', '#fe7e81']),
    'm embryonic': (['Pre-2C', '2C', '4C 8C', 'post-8C', 'Blastocyst', 
                     'Naive ESCs', 'Epiblast', 'ESCs', 'Late embryonic'], 
                    colors[3:12]),
}

# custom mpl.colormap for visualizing -1, 0, 1 descrete values
RdBu_bin = LinearSegmentedColormap.from_list('RdBu_binary', 
                                             [colors[18], '#f6f7f7', colors[14]], 
                                             3)

# dendrogram colors. By default set to all black.
dendrogram_colors = ['#000000', '#000000', '#000000', '#000000']

# dendrogram_colors = [colors[19], '#589909', '#0c9691', '#13bf63']
set_link_color_palette(dendrogram_colors[1:])

# default plot element sizes in inches
# HM_* args correspond to plot elements in both heatmap functions
HM_LEFT = 1.1                 # space left from the plot
HM_TOP = 1                    # space on top of the plot
HM_RIGHT = .2                 # space right from the plot
HM_BOTTOM = 1.6              # space at the bottom of the plot
HM_WSPACE = .04               # space between plot elements left and right
HM_HSPACE = .02               # space between plot elements top and bottom
HM_Y_COLORBAR = .04           # width of the colorbar on the left
HM_X_COLORBAR = .04           # width of the colorbar on the right
HM_DISTANCE_BAR = .08         # width of the distance bar on top of the heatmap
HM_Y_DENDROGRAM = .3          # size of the dendrogram on the y-axis
HM_X_DENDROGRAM = .3          # size of the dendrogram on the x-axis 
HM_SQUARE_SIZE = .1          # the size of one heatmap square

# G_HM_* args correspond to plot elements in the single-gene heatmap function
G_HM_SUMPLOT_SIZE = 1        # size of the summary plot
G_HM_UPDOWN_SPACE_SIZE = .4   # space between the up and down heatmaps

# CB_* args correspond to the colorbar dimensions in inches
CB_LEFT = .3                  # space left from the colorbar
CB_LEFT_SEC = 1.5             # space left from the second colorbar
CB_TOP = .35                  # space on top of the colorbar
CB_WIDTH = 1                  # width of the colorbar
CB_HEIGHT = .07               # height of the colorbar

# BP_* args correspond to plot elements in the barplot function
BP_LEFT = 1.6                 # space left from the bar plot 
BP_TOP = .7                   # space on top of the bar plot
BP_RIGHT = .5                 # space right from the bar plot
BP_BOTTOM = .5                # space on the bootom of the bar plot  
BP_Y_COLORBAR = .04           # width of the colorbar on the left
BP_BARSPACE = .8              # size of the bar plot (width)
BP_BARWIDTH_SIZE = .1         # width of the single bars

delta = r'$\Delta$ '
# plot labels for aggrevated information, per gene and specific distance bars 
EUCL_DIFF = ('{0}Expression similarity\n'
                    '[{0}Euclidean distance]'.format(delta))
# EUCL_DIFF_PROP = ('Proportion of changed expression similarity\n'
#                       '[proportional differential mean Euclidean distance]')
EUCL_NODIFF = ('Absolute expression similarity\n'
               '[Euclidean distance]')
# INTE_DIFF_NOPROP = ('Differential- & marker genes similarity\n'
#                         '[sum of matches(1) & mismatches(-1)]')
INTE_DIFF = ('Marker gene similarity\n'
             '[proportion of intersect]')

INTE_DIFF_SINGLE = ('Marker gene intersect\n'
                    '[matches and mismatches]')
INTE_DIFF_SINGLE_SUM = ('Proportion of intersect\n'
                      '[sum of matches(1) &\nmismatches(-1) /n]')

# EUCL_G_DIFF_NOPROP = ('Absolute expression similarity\n'
#                     '[differential Euclidean distance]')
# EUCL_G_DIFF_PROP = ('Proportion of changed expression similarity\n'
#                       '[proportional differential Euclidean distance]')
# EUCL_G_NODIFF = ('Absolute expression similarity\n'
#                '[absolute Euclidean distance]')
# INTE_G_DIFF_NOPROP = ('Differential genes similarity\n'
#                     '[target marker gene intersect]')

EUCL_DIST_BAR = ('Base Expression similarity\n'
                 '[Euclidean distance]')
INTE_DIST_BAR = ('Number of marker genes\n')
EUCL_G_DIST_BAR = ('Base expression similarity\n'
                   '[Euclidean distance]')


def _update_consts(constants):
    """take in the kwargs passed to the plot and override the config constans"""
    inv = [args for args in constants if args not in globals()]
    if inv:
        raise TypeError('{} are not recognized arguments.'.format(inv))
    globals().update(constants)