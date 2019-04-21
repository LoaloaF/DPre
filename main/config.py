from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from scipy.cluster.hierarchy import set_link_color_palette
import os


# DPI to from 100 to 300
mpl.rcParams['figure.dpi'] = 300.0
mpl.rcParams['savefig.dpi'] = 300.0
mpl.rcParams['font.size'] = FONTS = 2.66
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

DPRE_PATH =  os.path.dirname(__file__) + '/../'
SAVE_FORMAT = 'pdf'
UNDETECTED_MARKERGENES_BEHAVIOR =  'substitute'
UNDETECTED_MARKERGENES_BEHAVIOR =  'drop'
DESEQ2_P = .05
LOG_DEFAULT_TARGET_INIT_REDUCED = True

NORM_PROP_SINGGENE_EUCL_BEHAVIOUR = {
    'min_eucl_dist': .3,
    'treat': 'to_zero', #'drop'
    # 'treat': 'drop'
}

log_emphz = '=|=|=|=|=|=|=|=|====INITIATION====|=|=|=|=|=|=|=|=|='
log_plot = '=|=|=|=|=|=|=|=|====PLOT====|=|=|=|=|=|=|=|=|='


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
    'embryonic': 'fbcc8c', 
    'germ cells': '#f0f0a6',
    'neural crest': '#cap1cd', 
    'surface ectoderm': '#9edce4',
    'neuroectoderm': '#5888c2', 
    'mesoderm': '#64bb79', 
    'endoderm': '#843487', 
    'blood mesoderm': '#fe7e81', 
}

RdBu_bin = LinearSegmentedColormap.from_list(
            'RdBu_binary', [colors[18], '#f6f7f7', colors[14]], 3
)
# dendrogram_colors = [colors[19], '#589909', '#0c9691', '#13bf63']
dendrogram_colors = ['#000000', '#000000', '#000000', '#000000']
set_link_color_palette(dendrogram_colors[1:])



GSH_LEFT = .6
GSH_TOP = .6
GSH_RIGHT = .2
GSH_BOTTOM = .25
GSH_WSPACE = .04
GSH_HSPACE = .02

GSH_UP_DOWN_SPACE = .3
GSH_DRIVERS_COLORBAR = .03
GSH_GENES_COLORBAR = .03
GSH_DRIVERS_DENDROGRAM = .3
GSH_GENE_DENDROGRAM = .3
GSH_REQU_EFF_BAR = .04
GSH_SQUARESIZE = .006
GSH_SQUARESIZE = .05
GSH_SUMPLOT_SIZE = .9



TSH_LEFT = .6
TSH_TOP = .6
TSH_RIGHT = .01
TSH_BOTTOM = .7
TSH_WSPACE = .02
TSH_HSPACE = .02

TSH_DRIVERS_COLORBAR = .03
TSH_TARGETS_COLORBAR = .03
TSH_DRIVER_DENDROGRAM = .3
TSH_TARGET_DENDROGRAM = .3
TSH_REQU_EFF_BAR = .04
TSH_SQUARESIZE = .05



RSB_LEFT = .7
RSB_TOP = .3
RSB_RIGHT = .1
RSB_BOTTOM = .3

RSV_TARGETS_COLORBAR = .03
RSB_BARSPACE = 1
RSB_BARWIDTH_SIZE = .1
RSB_BETWEEN_BARS_SIZE = .04



CB_LEFT = .2
CB_LEFT_SEC = 1
CB_LEFT_TERT = 1.8
CB_TOP = .4
CB_WIDTH = .65
CB_HEIGHT = .04


AGG_EUCL_DIFF_NOPROP = ('mean change in expr. similarity\n'
                        '[differential mean eulc. dist.]')
AGG_EUCL_DIFF_PROP = ('prop. of changed expr. similarity\n'
                      '[prop. differential mean eucl. dist.]')
AGG_EUCL_NODIFF = ('mean abs. expr. similarity (line = base)\n'
                   '[mean abs. eucl. dist.]')
AGG_INTE_DIFF_NOPROP = ('differential genes similarity\n'
                        '[sum of matches(1) & mism.(-1)]')
AGG_INTE_DIFF_PROP = ('prop. of target markergene intersect\n'
                      '[sum of matches(1) & mism.(-1)]')




# def MAX_POSSIBLE_BAR_args(get_which):
#     # abs. mean euclidean distance --- number of intersect markergenes 
#     cmap =   ['afmhot', 'afmhot_r']
#     vmin =   [0,         100]
#     vmax =   [4,         1150]
#     aspect = ['auto',   'auto']

#     ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
#            for i in range(len(cmap))]
#     return ret[get_which] if get_which is not None else ret

# def TARGET_SIMILARITY_HEATMAP_args(get_which):
#     # euclidean distance - propotional euclidean distance
#     # markergene intersect - proportional markergene intersect
#     cmap =   ['RdBu_r', 'RdBu_r', 'RdBu_r', 'RdBu_r']
#     vmin =   [-.5,       -.15,     -65,      -.1]
#     vmax =   [.5,       .15,       65,       .1]
#     aspect = ['auto',  'auto',  'auto',    'auto']

#     ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
#            for i in range(len(cmap))]
#     return ret[get_which] if get_which is not None else ret


# def SUM_PLOT_args(get_which):
#     # margenetype 'up' --- margenetype 'down'
#     edgecolor = ['k',     'k']
#     height =    [.6,      .6]
#     xlim =      [(-.5, .5), (-2, 2)]
    
#     ret = [{'edgecolor': edgecolor[i], 'height': height[i], 'xlim': xlim[i]} 
#            for i in range(len(edgecolor))]
#     return ret[get_which] if get_which is not None else ret


# def SINGLE_GENE_HEATMAP_args(get_which):
#     # eucl. dist. --- eucl. dist. prop. --- eulcl. dist abs. --- intersect --- 
#     cmap =   ['RdBu_r', 'RdBu_r', 'afmhot_r', RdBu_bin]
#     vmin =   [-2.5,       -1,       0,          -1]
#     vmax =   [2.5,         1,       4,          1]
#     aspect = ['auto',   'auto',   'auto',     'auto']
    
#     ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
#            for i in range(len(cmap))]
#     return ret[get_which] if get_which is not None else ret
    

# def COLOR_BAR_args(get_which):
#     # y-colorbar --- x-colorbar
#     edgecolor = ['k', '']

#     ret = [{'edgecolor': edgecolor[i]} 
#            for i in range(len(edgecolor))]
#     return ret[get_which] if get_which is not None else ret



# def get_plot_args(of_plot, get_which=None):
#     if of_plot == 'max_possible_bar':
#         return MAX_POSSIBLE_BAR_args(get_which)
#     elif of_plot == 'accuracy_heatmap':
#         return ACCURACY_HEATMAP_args(get_which)
#     elif of_plot == 'target_similarity_heatmap':
#         return TARGET_SIMILARITY_HEATMAP_args(get_which)
#     elif of_plot == 'sum_plot':
#         return SUM_PLOT_args(get_which)
#     elif of_plot == 'single_gene_heatmap':
#         return SINGLE_GENE_HEATMAP_args(get_which)
#     elif of_plot == 'color_bar':
#         return COLOR_BAR_args(get_which)

