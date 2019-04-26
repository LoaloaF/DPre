from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from scipy.cluster.hierarchy import set_link_color_palette
import os

# DPI to from 100 to 300
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
DPRE_PATH =  os.path.dirname(__file__) + '/../'
SAVE_FORMAT = 'pdf'
UNDETECTED_MARKERGENES_BEHAVIOR = 'ignore'
UNDETECTED_MARKERGENES_BEHAVIOR = 'substitute'
DROP_TARGET_DETEC_THR = .15
DESEQ2_P = .05

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

default_targets_colors = {
    'embryonic': '#fbcc8c', 
    'germ cells': '#f0f0a6',
    'neural crest': '#c8b1ce', 
    'surface ectoderm': '#9edce4',
    'neuroectoderm': '#5888c2', 
    'mesoderm': '#64bb79', 
    'endoderm': '#843487', 
    'blood mesoderm': '#fe7e81', 
}

RdBu_bin = LinearSegmentedColormap.from_list('RdBu_binary', 
                                             [colors[18], '#f6f7f7', colors[14]], 
                                             3)
# dendrogram_colors = [colors[19], '#589909', '#0c9691', '#13bf63']
dendrogram_colors = ['#000000', '#000000', '#000000', '#000000']
set_link_color_palette(dendrogram_colors[1:])


HM_LEFT = .8
HM_TOP = .6
HM_RIGHT = .1
HM_BOTTOM = 1
HM_WSPACE = .04
HM_HSPACE = .02
HM_Y_COLORBAR = .04
HM_X_COLORBAR = .04
HM_REQU_EFF_BAR = .06
HM_Y_DENDROGRAM = .3
HM_X_DENDROGRAM = .3
HM_SQUARE_SIZE = .07

G_HM_SUMPLOT_SIZE = .9
G_HM_UPDOWN_SPACE_SIZE = .3

CB_LEFT = .2
CB_LEFT_SEC = 1
CB_TOP = .4
CB_WIDTH = .65
CB_HEIGHT = .06

BP_LEFT = 1
BP_TOP = .3
BP_RIGHT = .1
BP_BOTTOM = .3
BP_Y_COLORBAR = .04
BP_BARSPACE = .8
BP_BARWIDTH_SIZE = .07

AGG_EUCL_DIFF_NOPROP = ('mean change in expr. similarity\n'
                        '[differential mean eulc. dist.]')
AGG_EUCL_DIFF_PROP = ('prop. of changed expr. similarity\n'
                      '[prop. differential mean eucl. dist.]')
AGG_EUCL_NODIFF = ('mean abs. expr. similarity\n'
                   '[mean abs. eucl. dist.]')
AGG_INTE_DIFF_NOPROP = ('differential genes similarity\n'
                        '[sum of matches(1) & mism.(-1)]')
AGG_INTE_DIFF_PROP = ('prop. of target markergene intersect\n'
                      '[sum of matches(1) & mism.(-1)]')