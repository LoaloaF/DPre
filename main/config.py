from matplotlib.colors import LinearSegmentedColormap
import matplotlib as mpl
from scipy.cluster.hierarchy import set_link_color_palette


NORM_PROP_BAHVIOuR = {
    'min_abs_effect'
    'cap_type': 'proportion', #'hard_value'
    'cap': .05, #.3
    'change_prop_value': 'min_found', #'hard_cap', 'drop', 'none'
}
# NORM_PROP_BAHVIOuR = {
#     'cap_type': 'hard_value',
#     'cap': .4,
#     'change_prop_value': 0
# }

mpl.rcParams['font.size'] = 8
log_emphz = '=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|=|='

LOG_DEFAULT_TARGET_INIT_REDUCED = True

DPI = 300
SAVE_FORMAT = 'pdf'
# UNDETECTED_MARKERGENES_BEHAVIOR = 'drop'
UNDETECTED_MARKERGENES_BEHAVIOR = 'substitute'
DESEQ2_SIGNIFICANCE_THRESHOLD = .05

GREY1 = '#898989'
RED = '#8a0b25'
BLUE = '#0a3b70'

default_colors = ['#4daf4a', '#984ea3', '#ff7f00', '#ffff33', '#a65628', 
                  '#f781bf', '#999999', '#00ffff', '#5f0029', '#4f4e51',
                  '#b499ff', '#984ea3', '#00806c', '#27ff91', '#6a6c00', 
                  '#e1c78b', '#63c1fe', '#d90083', '#5a3500', '#42bba9',
                  '#b29a00','#f203ff','#004f74', '#ffa63a']

default_targets_colors = {
    'blood mesoderm': '#e1c78b', 
    'embryonic': '#63c1fe', 
    'endoderm': '#d90083', 
    'germ cells': '#5a3500',
    'mesoderm': '#42bba9', 
    'neural crest': '#b29a00', 
    'neuroectoderm': '#f203ff', 
    'surface ectoderm': '#004f74'  
}

RdBu_bin = LinearSegmentedColormap.from_list(
            'RdBu_binary', [RED, '#f6f7f7',BLUE], 3
)
dendrogram_colors = [GREY1, '#589909', '#0c9691', '#13bf63']
set_link_color_palette(dendrogram_colors[1:])



HMS_LEFT = .4
HMS_TOP = .6
HMS_RIGHT = .2
HMS_BOTTOM = .25
HMS_WSPACE = .04
HMS_HSPACE = .02

HMS_UP_DOWN_SPACE = .3
HMS_DRIVERS_COLORBAR = .03
HMS_GENES_COLORBAR = .03
HMS_DRIVERS_DENDROGRAM = .3
HMS_GENE_DENDROGRAM = .3
HMS_MAX_POSS_BAR = .04
HMS_SQUARESIZE = .006
HMS_SQUARESIZE = .05
HMS_SUMPLOT_SIZE = .9



HM_DRIVERS_COLORBAR = .03
HM_TARGETS_COLORBAR = .03
HM_DRIVER_DENDROGRAM = .3
HM_TARGET_DENDROGRAM = .3
HM_MAX_POSS_BAR = .04
HM_SQUARESIZE = .05

HM_LEFT = .6
HM_TOP = .6
HM_RIGHT = .01
HM_BOTTOM = .7
HM_WSPACE = .02
HM_HSPACE = .02



SE_YLABELSPACE = .9
SE_TARGETS_COLORBAR = .03
SE_BARSPACE = .8
SE_BARWIDTH_SIZE = .05
SE_BETWEEN_BARS_SIZE = .01
SE_XLIM = 3

SE_LEFT = .1
SE_TOP = .3
SE_RIGHT = .1
SE_BOTTOM = .2
SE_WSPACE = .02
SE_HSPACE = .2



CB_LEFT = .1
CB_LEFT_SEC = .9
CB_LEFT_TERT = 1.7
CB_TOP = .5
CB_WIDTH = .65
CB_HEIGHT = .06



def get_plot_args(of_plot, get_which=None):
    if of_plot == 'max_possible_bar':
        return MAX_POSSIBLE_BAR_args(get_which)
    elif of_plot == 'accuracy_heatmap':
        return ACCURACY_HEATMAP_args(get_which)
    elif of_plot == 'target_similarity_heatmap':
        return TARGET_SIMILARITY_HEATMAP_args(get_which)
    elif of_plot == 'sum_plot':
        return SUM_PLOT_args(get_which)
    elif of_plot == 'single_gene_heatmap':
        return SINGLE_GENE_HEATMAP_args(get_which)
    elif of_plot == 'color_bar':
        return COLOR_BAR_args(get_which)


def MAX_POSSIBLE_BAR_args(get_which):
    # abs. mean euclidean distance --- number of intersect markergenes 
    cmap =   ['afmhot', 'afmhot_r']
    vmin =   [0,         100]
    vmax =   [4,         1150]
    aspect = ['auto',   'auto']

    ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
           for i in range(len(cmap))]
    return ret[get_which] if get_which is not None else ret

def ACCURACY_HEATMAP_args(get_which):
    # euclidean distance - propotional euclidean distance
    # markergene intersect - proportional markergene intersect
    cmap =   ['viridis', 'viridis_r', 'viridis']
    vmin =   [0,         -1,           0]
    vmax =   [.8,         1,            65]
    aspect = ['auto',    'auto',       'auto']

    ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
           for i in range(len(cmap))]
    return ret[get_which] if get_which is not None else ret

def TARGET_SIMILARITY_HEATMAP_args(get_which):
    # euclidean distance - propotional euclidean distance
    # markergene intersect - proportional markergene intersect
    cmap =   ['RdBu_r', 'RdBu_r', 'RdBu_r', 'RdBu_r']
    vmin =   [-.5,       -.15,     -65,      -.1]
    vmax =   [.5,       .15,       65,       .1]
    aspect = ['auto',  'auto',  'auto',    'auto']

    ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
           for i in range(len(cmap))]
    return ret[get_which] if get_which is not None else ret


def SUM_PLOT_args(get_which):
    # margenetype 'up' --- margenetype 'down'
    edgecolor = ['k',     'k']
    height =    [.6,      .6]
    xlim =      [(-.5, .5), (-2, 2)]
    
    ret = [{'edgecolor': edgecolor[i], 'height': height[i], 'xlim': xlim[i]} 
           for i in range(len(edgecolor))]
    return ret[get_which] if get_which is not None else ret


def SINGLE_GENE_HEATMAP_args(get_which):
    # eucl. dist. --- eucl. dist. prop. --- eulcl. dist abs. --- intersect --- 
    cmap =   ['RdBu_r', 'RdBu_r', 'afmhot_r', RdBu_bin]
    vmin =   [-2.5,       -1,       0,          -1]
    vmax =   [2.5,         1,       4,          1]
    aspect = ['auto',   'auto',   'auto',     'auto']
    
    ret = [{'cmap': cmap[i], 'vmax': vmax[i], 'vmin': vmin[i], 'aspect': aspect[i]} 
           for i in range(len(cmap))]
    return ret[get_which] if get_which is not None else ret
    

def COLOR_BAR_args(get_which):
    # y-colorbar --- x-colorbar
    edgecolor = ['k', '']

    ret = [{'edgecolor': edgecolor[i]} 
           for i in range(len(edgecolor))]
    return ret[get_which] if get_which is not None else ret





def set_dpi(dpi, adjust_plot_sizes=True):
    dpi_rat = mpl.rcParams['figure.dpi'] /dpi
    mpl.rcParams['figure.dpi'] = dpi
    mpl.rcParams['savefig.dpi'] = dpi
    mpl.rcParams['font.size'] *=  dpi_rat
    mpl.rcParams['lines.linewidth'] *=  dpi_rat
    mpl.rcParams['hatch.linewidth'] *=  dpi_rat
    mpl.rcParams['patch.linewidth'] *=  dpi_rat
    mpl.rcParams['axes.linewidth'] *=  dpi_rat
    mpl.rcParams['xtick.major.size'] *= dpi_rat
    mpl.rcParams['xtick.major.width'] *= dpi_rat
    mpl.rcParams['ytick.major.size'] *= dpi_rat
    mpl.rcParams['ytick.major.width'] *= dpi_rat
    mpl.rcParams['grid.linewidth'] *= dpi_rat
    mpl.rcParams['xtick.major.pad'] *= dpi_rat
    mpl.rcParams['ytick.major.pad'] *= dpi_rat
set_dpi(DPI, False)

FS = mpl.rcParams['font.size']


