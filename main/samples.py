from DPre.main._differential import _differential
from DPre.main._logger import spacer, logger

class Samples(_differential):
"""class describing the data to explore similarity for

Samples can hold lists of diff. genes and expression data identifying a 
collection of samples. A control in the data extends functionality. 

Arguments:
    diff_genes (optional): directory with deseq2 output OR 2 directories with 
        up- and down genelists OR pandas.DataFrame. Defaults to None.
        For genelist input, filenames in the up- and down directory must be the
        same. When passing a DataFrame, the index should consist of all genes, 
        the columns of a pandas.MultiIndex with up or up and down at level 0 and 
        the element names at level 1. THe dtype is bool, diff. genes are marked 
        as True values. 
    expression (optional): directory with expression tsv file OR 
        pandas.DataFrame. Defaults to None.
        The .tsv input should have an ensg key in the first column or an ensg 
        column label. Columns `loc`, `name`, `tss_loc` and `strand` are removed
        automatically. The data should be exclusively numerical without NaN's.
        When passing a DataFrame, the data should be log2- and z-transformed. 
        The index should consist of all genes, the columns of a 
        pandas.MultiIndex with the element names at level 0, and `log2` & `z` at 
        level 1.
    ctrl (str, optional): The name of the control found in expression or diff. 
        genes. Defaults to None. Allowed to be missing in diff. genes.
    override_namematcher (bool, optional): When both diff_genes and expression 
        passed, this overrides the element names in diff. genes. Defaults to 
        False. When False, element names in diff. genes and expression are 
        expected to match perfectly.
    name (str, optional): Name label of the Samples. Defaults to 'Samples'. Used
        in logging and plot annotations.
    log: (bool, optional): Log the Samples initiation. Defaults to True.
    
Note:
    At least one of diff_genes and expression must be passed. When both are 
    passed, the inputs must have the same elment order.  

"""
    def __init__(self, diff_genes=None, expression=None, ctrl=None,
                 override_namematcher=False, name=None, log=True):
        self._ctrl = ctrl
        # call _differential __init__
        super().__init__(diff_genes=diff_genes, expression=expression, name=name,
                         override_namematcher=override_namematcher, log=log)
        
        # log the number of found differential genes
        if self._has_diff and log:
            spacer.info('')
            n = self._diff.sum().unstack(0).reindex(self._names).to_string()
            logger.info('Number of diff. genes: \n{}'.format(n))
        self._log_init(log, 'diff. genes')

    @property
    def _names_noctrl(self):
        ns = self._names
        if self._ctrl and self._ctrl in ns:
            ns.remove(self._ctrl)
        return ns

    def __repr__(self):
        return ('\n=|=|= {}-instance =|=|=\nname = {};\nelements = {};\n'
                'diff. genes data = {};\nexpression data = {} (n={});\n'
                .format(self._type_name, self._name, self._names, len(self),
                        self._has_diff, self._has_expr))