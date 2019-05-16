"""module storing package wide variables, colors and plot sizes. May be edited
by the user"""
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import set_link_color_palette

# change DPI from 100 to 300, adjust everything accordingly
mpl.rcParams['figure.dpi'] = 300.0
mpl.rcParams['savefig.dpi'] = 300.0
mpl.rcParams['font.size'] = FONTS = 3.5
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
# Target is dropped from the analysis. By default 15%.
DROP_TARGET_DETEC_THR = .15
# p value when dealing with deseq2 differential input
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
# labeled groups. For these preset_targets, colors are defined in colors.tsv.
preset_col_legend = {
    'm embryonic': (['pre-2C', '2C', '4C 8C', 'post-8C', 'Blastocyst', 
                     'Naive ESCs', 'Epiblast', 'ESCs', 'late embryonic'], 
                    colors[3:12]),
    'mouse': (['embryonic', 'germ cells', 'neural crest', 'surface ectoderm', 
               'neuroectoderm', 'mesoderm', 'endoderm', 'blood mesoderm'], 
              ['#fbcc8c', '#f0f0a6', '#c8b1ce', '#9edce4', '#5888c2',  
               '#64bb79',  '#843487',  '#fe7e81']),
    'human': (['embryonic', 'surface ectoderm', 'neuroectoderm', 'mesoderm',  
               'endoderm', 'blood mesoderm'],
              ['#fbcc8c', '#9edce4', '#5888c2', '#64bb79', '#843487', '#fe7e81']),
}

# custom mpl.colormap for visualizing -1, 0, 1 descrete values
RdBu_bin = LinearSegmentedColormap.from_list('RdBu_binary', 
                                             [colors[18], '#f6f7f7', colors[14]], 
                                             3)
# dendrogram colors. By default set to all black.
dendrogram_colors = ['#000000', '#000000', '#000000', '#000000']
# dendrogram_colors = [colors[19], '#589909', '#0c9691', '#13bf63']
set_link_color_palette(dendrogram_colors[1:])

# default plot elment sizes in inches
# HM_* args correspond to plot elements in both heatmap functions
HM_LEFT = 1.1
HM_TOP = 1
HM_RIGHT = .2
HM_BOTTOM = 1.35
HM_WSPACE = .04
HM_HSPACE = .02
HM_Y_COLORBAR = .04
HM_X_COLORBAR = .04
HM_DISTANCE_BAR = .06
HM_Y_DENDROGRAM = .3
HM_X_DENDROGRAM = .3
HM_SQUARE_SIZE = .07

# G_HM_* args correspond to plot elements in the single-gene heatmap function
G_HM_SUMPLOT_SIZE = .9
G_HM_UPDOWN_SPACE_SIZE = .4

# CB_* args correspond to the colorbar dimensions in inches
CB_LEFT = .3
CB_LEFT_SEC = 1.3
CB_TOP = .25
CB_WIDTH = .85
CB_HEIGHT = .06

# BP_* args correspond to plot elements in the barplot function
BP_LEFT = 1.2
BP_TOP = .7
BP_RIGHT = .5
BP_BOTTOM = .4
BP_Y_COLORBAR = .04
BP_BARSPACE = .8
BP_BARWIDTH_SIZE = .07

# Aggregated plot labels; used for all plots except single-gene annotation
AGG_EUCL_DIFF_NOPROP = ('Mean change in expr. similarity\n'
                        '[differential mean Eucl. dist.]')
AGG_EUCL_DIFF_PROP = ('Prop. of changed expr. similarity\n'
                      '[prop. differential mean Eucl. dist.]')
AGG_EUCL_NODIFF = ('Mean abs. expr. similarity\n'
                   '[mean abs. Eucl. dist.]')
AGG_INTE_DIFF_NOPROP = ('Diff.- & marker genes similarity\n'
                        '[sum of matches(1) & mism.(-1)]')
AGG_INTE_DIFF_PROP = ('Prop. of target markergene intersect\n'
                      '[sum of matches(1) & mism.(-1)]')

def _update_consts(constants):
    inv = [args for args in constants if args not in globals()]
    if inv:
        raise TypeError('{} are not recognized arguments.'.format(inv))
    globals().update(constants)