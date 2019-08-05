""" DPre is a bioinformatic tool that enables the user to explore cell type 
conversion/ differentiation experiments. DPre may take the expression or the 
differentially regulated genes of the samples to rate transcriptional similarity 
with a reference dataset, the targets. 
"""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
del os, sys

from .main.samples import samples
from .main.targets import targets
from .main._format_input import preset_targets
import main.config as config
from .main._dpre_util import plot_color_legend
from .main._dpre_util import annotate
from .main._dpre_util import get_ensgs
from .main._dpre_util import add_diff_genes_from_z