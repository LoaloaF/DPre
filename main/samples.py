from DPre.main._differential import _differential

class Samples(_differential):

    def __init__(self, diff_genes=None, expression=None, ctrl=None,
                 override_diff_names=False, name=None, log=True):
        super().__init__(diff_genes=diff_genes, expression=expression, ctrl=ctrl, 
                         override_diff_names=override_diff_names, name=name,
                         log=log)