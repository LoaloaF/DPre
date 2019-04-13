import os, sys
import pandas as pd
import glob

from DPre.main._logger import logger, spacer
import DPre.main.config as config
import DPre.main._dpre_util as util


def TARGET(get, sort=False, colors_from_file=False):
    spacer.info('\n\n')
    valid = ['all', 'blood mesoderm', 'embryonic', 'endoderm', 'germ cells',
            'mesoderm', 'neural crest', 'neuroectoderm', 'surface ectoderm']
    if get not in valid:
        logger.error('`{}` is not a valid default target. valid ones are {}'
                     .format(get, valid))
        sys.exit(1)

    get_dir = '{}default_targets/{}'.format(config.DPRE_PATH, get)
    try:         
        diff = pd.read_pickle('{}/differential.gzip'.format(get_dir))
        expr = pd.read_pickle('{}/expression.gzip'.format(get_dir))
    except FileNotFoundError as e:
        logger.error('Could not load data for default targets: {}'.format(e))
        sys.exit(1)

    if sort:
        diff.sort_index(axis=1, inplace=True)
        expr.sort_index(axis=1, inplace=True)
    
    from DPre.main.targets import Targets
    t = Targets(diff_genes=diff, expression=expr, name='{} lineages'.format(get),
                log=not config.LOG_DEFAULT_TARGET_INIT_REDUCED)
    if config.LOG_DEFAULT_TARGET_INIT_REDUCED:
        logger.info('Default target `{}` created, name: `{}`, elements: {}'
                    .format(get, t._name, len(t)))
    if colors_from_file:
        try:
            df_colors = pd.read_csv('{}/colors.tsv'.format(get_dir), sep='\t', 
                                    index_col=0)
            t.set_colors(dict(zip(df_colors.index, df_colors.color)), 
                         log=not config.LOG_DEFAULT_TARGET_INIT_REDUCED)
        except FileNotFoundError as e:
            logger.warning('No colors.tsv file found: {}'.format(e))
    elif get != 'all':
        # does not work for get == 'all'
        t.set_colors([config.default_targets_colors[get]],
                     log=not config.LOG_DEFAULT_TARGET_INIT_REDUCED)
    return t




def _validate_expression(expr, type_name, ctrl):
    isna = expr.isna()
    if isna.any().any():
        if not ((expr.columns.nlevels == 2) and not isna.xs('z', 1, 1).any().any()):
            spacer.error('\n')
            logger.error('Invalid expression data: data contains NaN values.')
            sys.exit(1)

    if 'object' in expr.dtypes.values:
        spacer.error('\n')
        logger.error('Invalid expression data: data contains elements of '
                     'datatype `object` (often text).')
        sys.exit(1)

    if type_name == 'Drivers' and not ctrl:
        spacer.error('\n')
        logger.error('For Drivers, the name of the contrl must be passed if '
                     'expression data is loaded.')
        sys.exit(1)
    
    elif ctrl and (ctrl not in expr.columns.unique(0)):
        spacer.error('\n')
        logger.error('The contrl name `{}` was not found in the passed '
                     'expression data.'.format(ctrl))
        sys.exit(1)
        
    if expr.columns.nlevels > 1:
        exp_idx = [(name, dt) for name in expr.columns.unique(0) 
                    for dt in ['log2', 'z']]
        idx = expr.columns.values.tolist()
        misma = list(filter(lambda i: i not in exp_idx, idx))
        if any(misma):
            spacer.error('')
            msg = ('\tInvalid expresion data. When passing data with '
                    'log2- and z-data, the columns must be a MultiIndex '
                    'in which level 0 holds the names: [`name1`, ...] '
                    'and level 1 the data types: [`log2`, `z`]. '
                    'Expected column indices ({}):\n\t\t{}\n\t Passed, '
                    'unexpected column indices ({}):\n\t\t{}'
                    .format(len(exp_idx), exp_idx, len(misma), misma))
            logger.error(msg)
            sys.exit(1)
        else:
            return expr
    else:
        return util._add_log2_z(expr)


def _format_diff_genes(diff_genes_dir, has_expr, genelists_mgtype='up'): 
    # check if path exists and contains .tsv files
    def check_path(direc):
        if not os.path.exists(direc):
            spacer.info('')
            logger.error('Could not change directory to {}/{}\nCheck the '
                        'path.'.format(os.getcwd(), direc))
            sys.exit(1)
        files = glob.glob(direc + '/*.tsv')
        if not files:
            spacer.info('')
            logger.error('No *.tsv files found in {}/{}\nCheck the path.'
                        .format(os.getcwd(), direc))
            sys.exit(1)

    # check if up and down genelist directories are compatible
    def check_up_down_genelists():
        if isinstance(diff_genes_dir, (list, tuple)):
            check_path(diff_genes_dir[0])
            check_path(diff_genes_dir[1])
            up_dir = glob.glob(diff_genes_dir[0]+'/*.tsv')
            down_dir = glob.glob(diff_genes_dir[1]+'/*.tsv')
            if len(up_dir) != len(down_dir):
                logger.error('Number of up- and down genelist files differ. '
                            'Found {} *.tsv files in up directory:\n{}/{}\n'
                            '{} *tsv files in down directory:\n{}/{}\n'
                            .format(len(up_dir), os.getcwd(), 
                                    diff_genes_dir[0], len(down_dir), 
                                    os.getcwd(), diff_genes_dir[1]))
                sys.exit(1)

            is_single = lambda n: (n not in up_dir) or (n not in down_dir)
            singles = list(filter(is_single, set((*up_dir, *down_dir))))
            if singles:
                logger.error('Names of up- and down genelist files differ. '
                            'Names only found in one of the two '
                            'directories ({}):\n{}'
                            .format(len(singles), singles))
                sys.exit(1)
            get_down_genelists = True
            return diff_genes_dir[0]
        else:
            get_down_genelists = False
            check_path(diff_genes_dir)
            return diff_genes_dir, get_down_genelists

    # check wether to process deseq2 files or genelist files
    def check_input_type():
        test_df = pd.read_csv(files[0], sep='\t')
        index_col = 'ensg'
        # else:
        if 'ensg' not in test_df.columns:
            index_col = test_df.columns[0]
            if not str(test_df[index_col][0]).startswith('ENS'):
                spacer.error('')
                logger.error('The *.tsv files holding the differential '
                            'keys do not have a column `ensg` nor do they '
                            'seem to have an ensg index in the first '
                            'column: {}, {}, ...'
                            .format(*test_df[index_col][:2].tolist()))
                sys.exit(1)
        deseq2_cols = ['Unnamed: 0', 'baseMean', 'log2FoldChange', 'lfcSE', 
                    'stat', 'pvalue', 'padj']
        inp_type =  'deseq2' if test_df.columns.tolist() == deseq2_cols \
                            else 'genelist ({})'.format(genelists_mgtype)
        return inp_type, index_col

    # glue files together to one DataFrame
    def merge_files():
        diffs = []
        for file in files:
            name = file[file.rfind('\\')+1:-4]
            usecols = (0, 2, 6) if inp_t == 'deseq2' else None
            d = pd.read_csv(file, sep='\t', index_col=index_col, usecols=usecols)
            if d.isna().any().any():
                spacer.info('')
                logger.warning('{} NaN values found and deleted in {}.tsv'
                            .format(d.isna().any(1).sum(), name))
            if inp_t == 'deseq2':
                sig = d.padj < config.DESEQ2_P
                up_reg = (sig & (d.log2FoldChange > 0)).rename(('up', name))
                down_reg = (sig & (d.log2FoldChange < 0)).rename(('down', name))
                diffs.append(pd.concat((up_reg, down_reg), axis=1))
            else:
                s = pd.Series(True, index=d.index, name=(genelists_mgtype, name))
                diffs.append(s)
        return diffs

    spacer.info('\n\n')
    direc, get_down_genelists = check_up_down_genelists()
    files = glob.glob(direc + '/*.tsv')
    inp_t, index_col = check_input_type()

    f_first = files[0][files[0].rfind('\\')+1:]
    f_last = files[-1][files[-1].rfind('\\')+1:]
    spacer.info('')
    logger.info('Formatting differential genes from {} files. {} *.tsv '
                'files in {}:\n{} ... {}'
                .format(inp_t, len(files), direc, f_first, f_last))
    if (inp_t == 'genelist (up)') and not has_expr:
        spacer.warning('')
        logger.warning('You are iniatiating the {} from genelist files ' 
                    'without passing expression data. Note that you '
                    'may therefore have a low value of detected genes.')

    diffs = merge_files()
    if inp_t == 'genelist (up)' and get_down_genelists:
        # recursion entry
        diffs.extend(_format_diff_genes(diff_genes_dir[1], has_expr, 'down'))
        # recursion over
    elif inp_t == 'genelist (down)':
        # recursion exit
        return diffs
    return pd.concat(diffs, axis=1, sort=True).fillna(False).sort_index(1)
