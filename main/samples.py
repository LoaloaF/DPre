from DPre.main._differential import _differential
from DPre.main._logger import spacer

class Samples(_differential):

    def __init__(self, diff_genes=None, expression=None, ctrl=None,
                 override_namematcher=False, name=None, log=True):

        self._ctrl = ctrl
        super().__init__(diff_genes=diff_genes, expression=expression, name=name,
                         override_namematcher=override_namematcher, 
                         log=log, diff_mg_logname='diff. genes')

        if log:
            spacer.info('\n\n')
            self._log_init()

    @property
    def _names_noctrl(self):
        ns = self._names
        if self._ctrl and self._ctrl in ns:
            ns.remove(self._ctrl)
        return ns

    def __repr__(self):
        return ('\n=|=|= {}-instance =|=|=\nname = {};\nelements = {};\n'
                'differential data = {};\nexpression data = {} (n={});\n'
                .format(self._type_name, self._name, self._names, len(self),
                        self._has_diff, self._has_expr))