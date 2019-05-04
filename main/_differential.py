import pandas as pd
import sys
import copy
import matplotlib.colors as mpl_colors

import DPre.main._dpre_util as util
import DPre.main._format_input as fmt_inp
import DPre.main.config as config
from DPre.main._logger import logger, spacer, log_init

class _differential:
    """base class of the core classes Targets and Drivers
    
    This class provides the basic functionality of an RNAseq dataset container 
    holding expression- and differential genelist data

    Attributes:
        _expr (pd.DataFrame): expression data (float). ENSG keys make up the,
            index, the columns have the element names at level 0 (MultiIndex); 
            each element has a log2 and z column.
        _diff (pd.DataFrame): differential gene data for Samples, markergene 
            data for Targets. Only the naming is differnt, the data structure
            for both is the same. The index is the same as the expression one,
            the columns are split in up and down markergenes, then the element
            names at level 1. A differentially regualted gene (or markergene
            respectively) has a True value.
        _has_diff (bool): instance was initiated with differential input
        _has_diff (bool): instance was initiated with expression input
        name (str): the name used in logging messages and plot annotations
        _colors (dict): a mapping of all elements name and their color
    """

    def __init__(self, diff_genes, expression, name, override_namematcher, log):
        # set the _expr attribute with experssion input if passed
        self._has_expr = False
        self._expr = expression
        if self._expr is not None:
            self._has_expr = True
            ctrl = self._ctrl if self._type_name == 'Samples' else None
            self._expr = fmt_inp._format_expr(expression, self._type_name, ctrl)

        # set the _diff attribute with diff_genes input if passed
        self._has_diff = False
        self._diff = diff_genes
        if self._diff is not None:
            self._has_diff = True
            if not isinstance(diff_genes, pd.DataFrame):
                self._diff = fmt_inp._format_diff_genes(diff_genes, 
                                                      type_name=self._type_name)
            if self._has_expr:
                self._diff = self._diff.reindex(self._expr.index).fillna(False)
                self._is_expr_diff_compatible(override_namematcher)
        
        # set basic attributes name and colors
        self.name = name if name else self._type_name
        self._colors = {}

    @property
    def names(self):
        """Get the element names
        
        Setting names:
            Pass a dictionary with keys of old names and values of new ones or
            a list of new names that should replace the old names in this order
            
        """
        name_from, at = (self._expr, 0) if self._has_expr else (self._diff, 1)
        return name_from.columns.unique(at).tolist()

    @names.setter
    def names(self, names):
        spacer.info('\n\n')
        prv_names = self.names
        if type(names) is not dict:
            if len(names) != len(prv_names):
                spacer.error('\n')
                logger.error('The passed list of element names ({}) must '
                             'have the same length as the current one ({}).'
                             .format(len(names), len(prv_names)))
                sys.exit(1)
            df = pd.DataFrame({'current names: ': prv_names, 
                               'passed names: ': names})
            logger.info('Ensure previous and passed names align:'
                        '\n{}\n'.format(df.to_string()))
            names = dict(zip(prv_names, names))
        else:
            inv = [k for k in names if k not in prv_names]
            if inv:
                spacer.error('\n')
                logger.error('Keys of the passed mapping are not containd in '
                             'current element names: {}'.format(inv))
                sys.exit(1)
        self._update_data_columns(names, func_name='rename')
        self._colors = {nn: self._colors[pn] for pn, nn in names.items() 
                        if pn in self._colors}
        if self._ctrl and (self._ctrl in names):
            self._ctrl = names[self._ctrl]

        logger.info('`{}` renamed. New names:\n{}'.format(self.name, self.names))
    
    @property
    def _detec_genes(self):
        """Get an pandas.Index of the detected genes in Targets/Samples"""
        genes_from = self._expr if self._has_expr else self._diff
        return genes_from.index

    @property
    def _type_name(self):
        """Get instance type name, `Targets` or `Samples`"""
        return self.__class__.__name__

    def __len__(self):
        return len(self.names)

    def set_name(self, name):
        """Set the Targets/Samples name"""
        self.name = name

    def get_colors(self, order=None):
        """Get the Targets/Samples colors of `order`

        Args:
            order (list, optional): A listing of element names in which the 
            element colors will be returned. Defaults to the instances element 
            order. 
        
        Note:
            Not contained element names raise an error. If the color for the 
            requested element hasn't been set, the color is set to white.
        
        """
        if not self._colors:
            self.set_colors(['#ffffff'], log=False)
            logger.warning('Colors have not been set for `{}`. All have been '
                           'set to default: white'.format(self.name))
        if order is not None:
            not_ctnd = list(filter(lambda o: o not in self.names, order)) 
            if not_ctnd:
                logger.error('Failed to get colors of `{}`. The passed order '
                             'contains elements other than the element names: '
                             '{}'.format(self.name, not_ctnd))
                sys.exit(1)
            not_set = list(filter(lambda o: o not in self._colors, order)) 
            if not_set:
                self.set_colors(dict.fromkeys(not_set, '#ffffff'), log=False)
                logger.warning('Failed to get some colors of `{}`. The passed '
                               'order contains elements without a set color: {}'
                               '. Colors for these were set to default: white.'
                               .format(self.name, not_set))
        else:
            if order is None:
                order = self.names
        return [self._colors[k] for k in order]

    def set_colors(self, colors, log=True):
        """Set the colors for elements in Targets/Samples

        Args:
            colors (dict:list): A mapping of element name to color or a list
                of colors that is assigned to the elements (in this order) until 
                its end. If the list has only one element, all elements are set 
                to this color. The color value must be a valid input for 
                matplotlib, i.e. 'w', '#ffffff',  (1, 1, 1) all refer to white.
            log (bool, optional): Log the set colors. Defaults to True.

        Note:
            Not contained element names are ignored. If the color is invalid,
            the element color is set to white.

        """
        if log:
            spacer.info('\n\n')
        if type(colors) is dict:
            inv = [col_key for col_key in colors if col_key not in self.names]
            if inv:
                logger.warning('Passed color mapping for `{}` contained keys '
                               'other than the element names: {}.\n'
                               .format(self.name, inv))
                colors = {k: colors[k] for k in colors if k in self.names}
        else:
            if len(colors) == 1:
                colors = dict.fromkeys(self.names, colors[0])
            else:
                colors = dict(zip(self.names, colors))

        inv_cs = set(tuple(filter(lambda c: not mpl_colors.is_color_like(c), 
                        list(colors.values()))))
        if any(inv_cs):
            [colors.update({k: 'w'}) for k, c in colors.items() if c in inv_cs]
            logger.warning('Invalid color values found: {}. Substituted with '
                           'default color: white.\n'.format(inv_cs))
        self._colors.update(colors)
        if log:
            str_cols = ',\n\t'.join([key+': '+str(col) 
                                    for key, col in self._colors.items()])
            logger.info('Colors set for `{}`:\n\t{}'.format(self.name, str_cols)) 
                    



    def reorder(self, order):
        """Reorder the elements in Targets/Samples inplace

            Args:
                order (list): the new order of elements

            Note:
                If not all current element names or new element names are 
                passed, an error is raised.
        """
        spacer.info('\n\n')
        not_ctnd = list(filter(lambda o: o not in self.names, order))
        missing = list(filter(lambda n: n not in order, self.names))
        if missing:
            spacer.error('')
            logger.error('The passed order misses current element names: {}'
                         .format(missing))
            sys.exit(1)
        if not_ctnd:
            spacer.error('')
            logger.error('Invalid element names. Passed elements not contained '
                         'in current element names:\n{}'.format(not_ctnd))
            sys.exit(1)
        self._update_data_columns(order)
        logger.info('`{}` reordered. New order of elements:\n{}'
                    .format(self.name, self.names))

    def slice_elements(self, elements, name=None, log=True):
        """Slice the Targets/Samples to a specific list of elements. Return a 
            copy of the original.

        Args:
            elements (list): the list of elements to slice. 
        
        Returns:
            sliced: the sliced Targets/Samples instance
        
        Note:
            If the passed names are not found in the current element names, 
            an error is raised.
        """
        if log:
            spacer.info('\n\n')
        if not elements or not len(elements):
            spacer.error('')
            logger.error('The list of elements cannot be empty.')
            sys.exit(1)
        not_ctnd = list(filter(lambda e: e not in self.names, elements))
        if not_ctnd:
            spacer.error('')
            logger.error('Invalid element names. Passed elements not contained '
                         'in current element names:\n{}'.format(not_ctnd))
            sys.exit(1)
        if self._type_name == 'Samples': 
            if self._ctrl and (self._ctrl not in elements):
                self._ctrl = None

        sliced = copy.copy(self)
        sliced._update_data_columns(elements)
        [sliced._colors.pop(k, None) for k in sliced.names if k not in elements]
        sliced.name = name if name else sliced.name

        if log:
            logger.info('`{}` sliced:'.format(self.name))
            spacer.info('')
            diff_n = 'markergenes' if self._type_name == 'Targets' else \
                     'diff. genes' 
            sliced._log_init(log, diff_n)
        return sliced

    def _is_expr_diff_compatible(self, override_namematcher=False):
        """Check if expression and differential input share element labels

        Args:
            override_namematcher: Whether a misalignment between expression- and
                differnetial names should be ignored. Useful when names refer to
                the same data but labels differ partially.
        """
        expr_ns = self.names
        diff_ns = self._diff.columns.unique(-1).tolist()
        diff_n = 'markergenes' if self._type_name == 'Targets' else 'diff. genes'
        spacer.info('\n\n')
        if len(expr_ns) != len(diff_ns):
            if (len(expr_ns) == len(diff_ns)+1) and self._type_name == 'Samples':
                msg = ('{} ({}) has one element less than expression ({}). '
                       .format(diff_n, len(diff_ns), len(expr_ns)))
                if self._ctrl:
                    # assume it is the contrl missing
                    msg += ('An empty element `{}` (control) will be added to '
                            'diff. genes.'.format(self._ctrl))
                    logger.info(msg)

                    # add control to match with expression element names
                    self._diff[('up', self._ctrl)] = False
                    self._diff[('down', self._ctrl)] = False
                    self._diff.sort_index(axis=1, inplace=True)
                    diff_ns = self._diff.columns.unique(-1).tolist()
                else:
                    msg += ('If the expression data has a control that is '
                            'missing in diff. genes, you can resolve this by '
                            'passing the control name for Samples initiation.')
                    logger.info(msg)
                    sys.exit(1)

            else:
                logger.error('Passed expression ({}) and {} ({}) do not have '
                             'the same number of elements. Check input.'
                             .format(len(expr_ns), diff_n, len(diff_ns)))
                sys.exit(1)
        
        align = [e_n == d_n for e_n, d_n in zip(expr_ns, diff_ns)]
        df = pd.DataFrame({'expression names': expr_ns, 
                           diff_n + ' names': diff_ns, 
                           'match': align})
        if all(align) or override_namematcher:
            spacer.info('')
            msg = diff_n + ' names have been overriden by expression names. '
            lvls = self._diff.columns.unique(0), expr_ns
            self._diff.columns = pd.MultiIndex.from_product(lvls)
            if override_namematcher:
                logger.warning('{}CAUTION: {}- and expression names '
                               'do not match for all elements:\n{}\nMake sure '
                               'data aligns to avaid mislabeling!'
                               .format(msg, diff_n, df.to_string()))
        else:
            spacer.error('')
            logger.error(('{0}- and expression element names do '
                          'not match:\n{1}\nRename elements in expression '
                          '/{0} input files or override the '
                          '{0} names by setting `override_namematcher` '
                          'to True.'.format(diff_n, df.to_string())))
            sys.exit(1)

    def _log_init(self, log):
        """Check if expression or differential was passed, then log the 
        initiated Targets/Samples result

        Args:
            log: log initiation, otherwise only check input
            diff_n: differential name, 'markergenes' for Targets diff. genes'
                for Samples
        """
        if not self._has_diff and not self._has_expr:
            spacer.error('')
            cmd_pref = self._type_name.lower()+'_' if log == 'from_cmd' else ''
            logger.error('Failed to init {}:\nAt least one of `{}expression` '
                         'and `{}` must be passed.'
                         .format(self._type_name, cmd_pref, diff_n))
            sys.exit(1)
        if not log:
            return
        # assemble log message by checking the various inputs
        diff_n = 'markergenes' if self._type_name == 'Targets' else 'diff. genes'
        if self._has_diff:
            d = self._diff.any(axis=1).sum()
            n_diff_descr = (', of which {} {}'.format(d, diff_n))
        elif self._has_expr:
            n_diff_descr = ''

        if self._type_name == 'Targets':
            diffmgs = 'Markergenes'
            trg_smp_arg = 'Down markergenes loaded: {}'.format(self._down_mgs)
        elif self._type_name == 'Samples':
            diffmgs = 'Differential genes'
            trg_smp_arg = ('Control passed: {}'
                           .format(self._ctrl if self._ctrl else False))
        msg = ('{}\nNew {}-instance created: `{}`\n\t{} ({}):\n\t\t{}\n\t'
               'Detected genes: {}{}\n\t{} loaded: {}\n\t'
               'Expression data loaded: {}\n\t{}'
               .format(log_init, self._type_name, self.name,
                       self._type_name, len(self), ',\n\t\t'.join(self.names), 
                       len(self._detec_genes), n_diff_descr, diffmgs, 
                       self._has_diff, self._has_expr, trg_smp_arg))
        logger.info(msg)


    def _update_data_columns(self, names, func_name='reindex'):
        """Reindexer and renamer of the main datatypes
        
        Args:
            names: a mapper (dict) for func_name 'rename' or a list of contained
                names for reindex
            func_name: 'reindex' or 'rename'. The operation to execute.
        """
        if func_name == 'reindex':
            args = {'labels': names, 'axis': 1,}
        else:
            args = {'mapper': names, 'axis': 1,}

        if self._has_expr:
            if func_name == 'reindex':
                self._expr = self._expr.reindex(**args, level=0)
            elif func_name == 'rename':
                self._expr = self._expr.rename(**args, level=0)
                
            if self._type_name == 'Targets':
                if func_name == 'reindex':
                    self._expr_mgs = self._expr_mgs.reindex(**args, level=1)
                    self._expr_mgs.dropna(how='all', inplace=True)
                elif func_name == 'rename':
                    self._expr_mgs = self._expr_mgs.rename(**args, level=1)
                
        if self._has_diff:
            if func_name == 'reindex':
                self._diff = self._diff.reindex(**args, level=1)
            elif func_name == 'rename':
                self._diff = self._diff.rename(**args, level=1)