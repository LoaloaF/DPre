import sys, os

from DPre.main._differential import _differential
import DPre.main._dpre_util as util


                                                                       

class Drivers(_differential):

    n_insts = 0
    def __init__(self, diff_genes=None, expression=None, ctrl=None,
                 override_diff_names=False, name=None, log=True):

        super().__init__(diff_genes=diff_genes, expression=expression, ctrl=ctrl, 
                         override_diff_names=override_diff_names, name=name,
                         log=log)
        
        if self._has_diff:
            self._diff_eff = util._diff_to_int_updown_notation(self._diff)

        if self._has_expr:
            make_diff = lambda drv: drv - self._expr[self._ctrl].values
            expr_eff = self._expr.groupby(level=0, axis=1).apply(make_diff)
            self._expr_eff = util._add_mg_types(expr_eff, down=True)

        self.__class__.n_insts += 1
    
