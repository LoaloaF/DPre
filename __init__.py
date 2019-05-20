""" Initiation of DPre and import checking. Runs when `import DPre` is called. 
Defines which modules and functions are accessable through DPre.*
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