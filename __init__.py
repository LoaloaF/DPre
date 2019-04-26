import os, sys
sys.path.insert(0, os.path.abspath(os.path.dirname(__file__)))
del os, sys

from .main.samples import Samples
from .main.targets import Targets
from .main._format_input import TARGET
import main.config as config
from .main._dpre_util import color_legend