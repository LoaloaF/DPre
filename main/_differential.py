import pandas as pd
import sys
import copy
import matplotlib.colors as mpl_colors
import matplotlib.pyplot as plt

import DPre.main._dpre_util as util
import DPre.main._format_input as _format_input
import DPre.main.config as config
from DPre.main._logger import logger, spacer, log_init

class _differential:
    def __init__(self, diff_genes,  expression, ctrl, override_diff_names,
                 name, log):

        self._name = name if name else self._type_name
        self._ctrl = ctrl
        self._colors = {}
        self._has_expr = False
        self._has_diff = False

        self._expr = expression
        if self._expr is not None:
            self._has_expr = True
            self._expr = _format_input._validate_expression(expression, 
                                                    self._type_name, self._ctrl)
            self._min_zval = self._expr.xs('z', level=1, axis=1).min().min()
            self._max_zval = self._expr.xs('z', level=1, axis=1).max().max()

        self._diff = diff_genes
        if self._diff is not None:
            self._has_diff = True
            if not isinstance(diff_genes, pd.DataFrame):
                self._diff = _format_input._format_diff_genes(diff_genes, 
                                                              self._has_expr)

            if self._type_name == 'Targets':
                if not self._down_mgs:
                    self._diff = self._diff.xs('up', 1, 0, False)
                if 'down' not in self._diff.columns.unique(0):
                    self._down_mgs = False
            if log:
                spacer.info('')
                logger.info('Number of differential genes: \n{}'
                            .format(self._diff.sum().unstack(level=0).to_string()))
                        
            if self._has_expr:
                self._diff = self._diff.reindex(self._expr.index).fillna(False)
                self._is_expr_diff_compatible(override_diff_names, log)
        elif self._type_name == 'Targets':
            self._down_mgs = False
            logger.warning('The targets `{}` are initiated without marker '
                           'genes. Note that comparing against all genes can '
                           'lead to low accuracy for defining transcr. '
                           'similarity.'.format(self._name))

        
        if not self._has_expr and not self._has_diff:
            spacer.error('')
            logger.error('One of `expression` and `diff_genes` must be passed.')
            sys.exit(1)
        if log:
            spacer.info('\n\n')
            self._log_init()

    @property
    def _names_noctrl(self):
        ns = self._names
        if self._ctrl and self._ctrl in ns:
            ns.remove(self._ctrl)
        return ns
    
    @property
    def _names(self):
        name_from, at = (self._expr, 0) if self._has_expr else (self._diff, 1)
        return name_from.columns.unique(at).tolist()

    @property
    def _detec_genes(self):
        genes_from = self._expr if self._has_expr else self._diff
        return genes_from.index

    @property
    def _type_name(self):
        return self.__class__.__name__

    def __len__(self):
        return len(self._names)
    
    def __repr__(self):
        return ('\n=|=|= {}-instance =|=|=\nname = {};\nelements = {};\n'
                'differential data = {};\nexpression data = {} (n={});\n'
                .format(self._type_name, self._name, self._names, len(self),
                        self._has_diff, self._has_expr))



    def set_name(self, name):
        self._name = name

    def get_name(self, name):
        return self._name

    def get_names(self):
        return self._names

    def set_names(self, names):
        spacer.info('\n\n')
        prv_names = self._names
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
        
        self._update_data_columns(names, func_name='rename')
        self._colors = {nn: self._colors[pn] for pn, nn in names.items() 
                        if pn in self._colors}
        if self._ctrl and (self._ctrl in names):
            self._ctrl = names[self._ctrl]

        logger.info('`{}` renamed. New names:\n{}'.format(self._name, self._names))

        

    

    def get_colors(self, order=None):
        if not self._colors:
            self.set_colors(['#ffffff'], log=False)
            logger.warning('Colors have not been set for `{}`. All have been '
                           'set to default: white'.format(self._name))
        if order is not None:
            not_ctnd = list(filter(lambda o: o not in self._names, order)) 
            if not_ctnd:
                logger.error('Failed to get colors of `{}`. The passed order '
                             'contains elements other than the element names: '
                             '{}'.format(self._name, not_ctnd))
                sys.exit(1)
            not_set = list(filter(lambda o: o not in self._colors, order)) 
            if not_set:
                self.set_colors(dict.fromkeys(not_set, '#ffffff'), log=False)
                logger.warning('Failed to get some colors of `{}`. The passed '
                               'order contains elements without a set color: {}. '
                               'Colors for these were set to default: white.'
                               .format(self._name, not_set))
        else:
            if order is None:
                order = self._names
        return [self._colors[k] for k in order]

    def set_colors(self, colors, log=True):
        if log:
            spacer.info('\n\n')
        if type(colors) is dict:
            not_ctnd = list(filter(lambda k: k not in self._names, 
                                   list(colors.keys())))
            if not_ctnd:
                logger.warning('Passed color mapping for `{}` contained keys '
                               'other than the element names: {}.\n'
                               .format(self._name, not_ctnd))
                colors = {k: colors[k] for k in colors if k in self._names}
        else:
            if len(colors) == 1:
                colors = dict.fromkeys(self._names, colors[0])
            else:
                colors = dict(zip(self._names, colors))

        inv_cs = set(tuple(filter(lambda c: not mpl_colors.is_color_like(c), 
                        list(colors.values()))))
        if any(inv_cs):
            [colors.update({k: '#ffffff'}) for k, c in colors.items() if c in inv_cs]
            logger.warning('Invalid color values found: {}. Substituted with '
                           'default color: white.\n'.format(inv_cs))
        self._colors.update(colors)
        if log:
            str_cols = ',\n\t'.join([key+': '+str(col) 
                                    for key, col in self._colors.items()])
            logger.info('Colors set for `{}`:\n\t{}'.format(self._name, str_cols)) 
                    



    def reorder(self, order):
        spacer.info('\n\n')
        not_ctnd = list(filter(lambda o: o not in self._names, order))
        missing = list(filter(lambda n: n not in order, self._names))
        if missing:
            spacer.error('')
            logger.error('The passed order misses current element names: {}'
                         .format(missing))
        if not_ctnd:
            spacer.error('')
            logger.error('Invalid element names. Passed elements not contained '
                         'in current element names:\n{}'.format(not_ctnd))
        if missing or not_ctnd:
            sys.exit(1)
        self._update_data_columns(order)
        logger.info('`{}` reorderd. New order of elements:\n{}'
                    .format(self._name, self._names))

    def slice_elements(self, elements, name=None, log=True):
        if log:
            spacer.info('\n\n')
        if not len(elements):
            spacer.error('')
            logger.error('The list of elements cannot be empty.')
            sys.exit(1)
        not_ctnd = list(filter(lambda e: e not in self._names, elements))
        if not_ctnd:
            spacer.error('')
            logger.error('Invalid element names. Passed elements not contained '
                         'in current element names:\n{}'.format(not_ctnd))
            sys.exit(1)
        if self._ctrl and (self._ctrl not in elements):
            if (self._type_name == 'Samples') and self._has_expr:
                spacer.warning('')
                logger.warning('When slicing a Samples-instance holding '
                               'expression data, the control must be kept. `{}`'
                               ' was added to passed list of `elements`.'
                               .format(self._ctrl))
                elements.append(self._ctrl)
            else:
                self._ctrl = None

        slc = copy.copy(self)
        slc._update_data_columns(elements)
        [slc._colors.pop(k, None) for k in slc._names if k not in elements]
        slc._name = name if name else slc._name

        if log:
            logger.info('`{}` sliced:'.format(self._name))
            spacer.info('')
            slc._log_init()
        return slc

















    

    


    def _is_expr_diff_compatible(self, override_diff_names=False, log=True):
        expr_names = self._names
        diff_names = self._diff.columns.unique(-1).tolist()
        if log:
            spacer.info('\n\n')
        
        if len(expr_names) != len(diff_names):
            if (len(expr_names) == len(diff_names) +1):
                msg = ('Differential ({}) has one element less than expression '
                       '({}). '.format(len(diff_names), len(expr_names)))
                if self._ctrl:
                    msg += ('An empty element `{}` (control) will be added to '
                            'differential.'.format(self._ctrl))
                    logger.info(msg)

                    # add control to match with expression element names
                    self._diff[('up', self._ctrl)] = False
                    if 'down' in self._diff.columns.unique(0):
                        self._diff[('down', self._ctrl)] = False
                    self._diff.sort_index(axis=1, inplace=True)
                    diff_names = self._diff.columns.unique(-1).tolist()

                else:
                    msg += ('If the control is missing in the differential '
                            'elements, you can pass its name to append it '
                            'automatically.')
                    logger.info(msg)
                    sys.exit(1)

            else:
                logger.error('Passed expression ({}) and differential ({}) do '
                             'not have the same number of elements. Check input.'
                             .format(len(expr_names), len(diff_names)))
                sys.exit(1)
        
        align = [e_n in d_n for e_n, d_n in zip(expr_names, diff_names)]
        df = pd.DataFrame({'expression names': expr_names, 
                            'differential names': diff_names, 
                            'expr. substring of diff.': align})
        if all(align) or override_diff_names:
            spacer.info('')
            msg = 'Differential names have been overriden by expression names. '
            lvls = self._diff.columns.unique(0), expr_names
            self._diff.columns = pd.MultiIndex.from_product(lvls)
            if all(align):
                if log:
                    logger.info('Expression and differential names match: \n{}'
                                '\n{}'.format(df.to_string(), msg))
            elif override_diff_names:
                logger.warning('{}CAUTION: Differential- and expression names '
                               'do not match for all elements:\n{}\nMake sure '
                               'data aligns to avaid mislabeling!'
                               .format(msg, df.to_string()))
        else:
            spacer.error('')
            logger.error(('Differential- and expression element names do '
                          'not match:\n{}\nRename elements in expression '
                          '/differential input files or override the '
                          'differential names by setting `override_diff_names` '
                          'to True.'.format(df.to_string())))
            sys.exit(1)


    def _log_init(self):
        if self._has_diff:
            ud = self._diff.any(axis=1).sum()
            uniq_diff = (', of which {} unique differential'.format(ud))
        else:
            uniq_diff = ''

        msg = ('{}\nNew {}-instance created: `{}`\n\t{} ({}):\n\t\t{}\n\t'
               'Detected genes: {}{}\n\tDifferential genes loaded: {}\n\t'
               'Expression data loaded: {}'
               .format(log_init, self._type_name, self._name,
                       self._type_name, len(self), ',\n\t\t'.join(self._names), 
                       len(self._detec_genes), uniq_diff, self._has_diff, 
                       self._has_expr))
        if self._type_name == 'Samples':
            msg += '\n\tPassed control: {}'.format(bool(self._ctrl))
        logger.info(msg)


    def _update_data_columns(self, names, func_name='reindex'):
        args = {'labels': names, 'axis':1,}
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