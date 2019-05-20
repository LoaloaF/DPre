from DPre.main._differential import _differential
from DPre.main._logger import spacer, logger

class samples(_differential):
    """The data to explore similarity for.

    samples can hold lists of differential genes and expression data identifying 
    a collection of experimental samples. A control in the data enables the 
    computation of differential transcriptional similarity.

    Arguments:
        diff_genes (optional): Directory with deseq2 output, 2 directories
            with up- and down genelists or pandas.DataFrame. 
            For gene list input, filenames in the up- and down directory must be 
            the same. When passing a DataFrame, the index should consist of ensg
            keys, the columns of a pandas.MultiIndex with 'up' & 'down' at 
            level 0 and the element names at level 1. The dtype is bool, diff. 
            genes are marked as True values. 
        expression (optional): directory with expression TSV file or 
            pandas.DataFrame. Defaults to None.
            The TSV input should have an ensg key in the first column or an ensg 
            column label. Columns `loc`, `name`, `tss_loc` and `strand` are 
            removed automatically. The data should be exclusively numerical 
            without NaN's. When passing a DataFrame, the data can be log2- 
            and z-transformed with an ensg key index and pandas.MultiIndex 
            columns with the element names at level 0, and `log2` & `z` at 
            level 1. 
        ctrl (str, optional): The name of the control found in 'expression'.
            Defaults to None. Allowed to be missing in 'diff_genes'.
        override_namematcher (bool, optional): When both 'diff_genes' and 
            'expression' passed, this overrides the element names in diff_genes. 
            Defaults to False. When False, element names in 'diff_genes' and 
            'expression' are expected to match perfectly.
        name (str, optional): Name label of the samples. Defaults to 'Samples'. 
            Used in logging and plot annotations.
        log: (bool, optional): Log the samples initiation. Defaults to True.
        
    Note:
        At least one of 'diff_genes' and 'expression' must be passed. When both 
        are passed, the inputs must have the same element order. Gene list
        data is automatically alphabetically sorted, hence the expression order
        should concur with this.
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
            n = self._diff.sum().unstack(0).reindex(self.names).to_string()
            logger.info('Number of diff. genes: \n{}'.format(n))
        if log:
            spacer.info('\n\n')
        self._log_init(log)

    @property
    def _names_noctrl(self):
        """Return the element names without the control"""
        ns = self.names
        if self._ctrl and self._ctrl in ns:
            ns.remove(self._ctrl)
        return ns

    def __repr__(self):
        """Get a readable summary of the samples instance"""
        return ('\n=|=|= samples-instance =|=|=\nname = {};\nelements = {};\n'
                'n = {};\ndiff. genes data = {};\nexpression data = {};\n'
                .format(self.name, self.names, len(self), self._has_diff, 
                        self._has_expr))