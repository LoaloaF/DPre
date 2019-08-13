"""Handles user data input and initiation of preset targets"""
import os, sys
import glob
import pandas as pd

from DPre.main._logger import logger, spacer
import DPre.main.config as config
import DPre.main._dpre_util as util

def preset_targets(get, sort=False, preset_colors=True, 
                   color_legend_filename=True, color_legend_ncols=1):
    """Generate one of the predefined targets instances and return it. 
    
        Pick a reference dataset for comparison. Mouse (Hutchins et al. 2017, 
        NAR) and Human (Abugessaisa et al. 2017, FANTOM5 project) are included. 
        Specific doamins can be picked for both species references. If the 
        targets are initiated with 'preset_colors', a color legend is generated 
        and saved in the current working directory. Custom presets can be 
        created by adding a folder (with an 'm' or 'h' prefix) to 
        DPre/preset_targets. 

    Args:
        get (str): the desired preset. Valid options are 'mouse', 'human', 
            'm embryonic', 'm germ cells', 'm neural crest', 
            'm surface ectoderm', 'm neuroectoderm', 'm mesoderm', 'm endoderm', 
            'm blood mesoderm', 'h surface ectoderm', 'h neuroectoderm', 
            'h mesoderm', 'h endoderm', 'h blood mesoderm'; m = mouse, h = human
        sort (bool, optional): Sort the loaded element names alphabetically. 
            Defaults to False.
        preset_colors (bool, optional): Tries to initiate the targets with preset 
            colors either from colors.tsv in the respective preset directory or 
            when not found from config.preset_targets_colors. Defaults to True.
        color_legend_filename (bool, str, optional): The filename when a preset
            color legend is drawn from config.preset_col_legend. When True, a
            filename is inferred from the targets name and config.SAVE_FORMAT, 
            a str is set as the filename. Defaults to True. When None, the color
            legend is not drawn.
        color_legend_ncols (int, optional): Number of columns in the color 
            legend. Defaults to 1.
    
    Returns:
        t: the preset targets instance
    """
    path = os.path.dirname(__file__)
    # any folder in DPre/preset_targets is potentially valid
    valid = os.listdir(os.path.join(path, '..', 'preset_targets'))
    if get not in valid:
        spacer.info('')
        logger.error('`{}` is not a valid preset target. Valid ones are {}'
                     .format(get, valid))
        sys.exit(1)
    
    # try to get .gzip markergene and expression input, if not found try .tsv
    get_dir = '{}/../preset_targets/{}'.format(path, get)
    expr = mgs = None
    if os.path.exists('{}/markergenes.gzip'.format(get_dir)):
        mgs = pd.read_pickle('{}/markergenes.gzip'.format(get_dir))
    elif os.path.exists('{}/markergenes.tsv'.format(get_dir)):
        mgs = pd.read_csv('{}/markergenes.tsv'.format(get_dir), sep='\t', 
                          index_col=0, header=[0,1])
        mgs.to_pickle('{}/markergenes.gzip'.format(get_dir))
    
    if os.path.exists('{}/expression.gzip'.format(get_dir)):
        expr = pd.read_pickle('{}/expression.gzip'.format(get_dir))
    elif os.path.exists('{}/expression.tsv'.format(get_dir)):
        expr = pd.read_csv('{}/expression.tsv'.format(get_dir), sep='\t', 
                           index_col=0, header=[0,1])
        expr.to_pickle('{}/expression.gzip'.format(get_dir))
    elif os.path.exists('{}/expression_h1.gzip'.format(get_dir)):
        expr1 = pd.read_pickle('{}/expression_h1.gzip'.format(get_dir))
        expr2 = pd.read_pickle('{}/expression_h2.gzip'.format(get_dir))
        expr = pd.concat([expr1, expr2], axis=1)
        expr.to_pickle('{}/expression.gzip'.format(get_dir))

    if sort:
        mgs.sort_index(axis=1, inplace=True)
        expr.sort_index(axis=1, inplace=True)
    
    # explicit part of the script that might need adjustmet for cumstom presets
    if get == 'human':
        args = {'name': 'human FANTOM5 library', 'species': 'human'}
    elif get == 'mouse':
        args = {'name': 'mouse lineages', 'species': 'mouse'}
    elif get.startswith('h '):
        args = {'name': get[2:]+' lineage', 'species': 'human'}
    elif get.startswith('m '):
        args = {'name': get[2:]+' lineage', 'species': 'mouse'}
    # init targets
    from DPre.main.targets import targets
    t = targets(markergenes=mgs, expression=expr, log=False, **args)
    logger.info('Default targets `{}` created, name: `{}`, elements: {}'
                .format(get, t.name, len(t)))

    # try to get colors first through a file, then through config
    if preset_colors:
        try:
            df_colors = pd.read_csv('{}/colors.tsv'.format(get_dir), sep='\t', 
                                    index_col=0)
            t.set_colors(dict(zip(df_colors.index, df_colors.color)), log=False)
        except FileNotFoundError:
            if get in config.preset_targets_colors:
                t.set_colors([config.preset_targets_colors[get]], log=False)
            else:
                logger.warning('No colors found for preset targets {}'
                               .format(get))
        # draw a colorlegend if defined in config
        if get in config.preset_col_legend and color_legend_filename:
            filename = get+'_color_legend' if color_legend_filename == True \
                       else color_legend_filename
            util.plot_color_legend(*config.preset_col_legend[get], 
                                   ncolumns=color_legend_ncols, filename=filename)
    return t

def _format_expr(expr, type_name, ctrl=None):
    """ Take user expression input validate and format

    If a TSV file is passed, read the expresion file as a DataFrame. Check
    if the DataFrame has a valid format. If the control is passed, check if it's
    found in expression. Finally, generate and add the log2- and z-transformed 
    data.

    Args:
        expr: Filename or Dataframe. The data to check.
        type_name: 'targets' or 'samples', depending on caller
        ctrl: Control name, only passed when called from samples
    
    Returns:
        expr: Expression DataFrame with log2- and z-transformed data at column 
            level 1

    """
    if not isinstance(expr, pd.DataFrame):
        if not os.path.exists(expr):
            spacer.info('')
            logger.error('Invalid path: {}\n'.format(os.path.abspath(expr)))
            sys.exit(1)

        expr = pd.read_csv(expr, sep='\t')
        if 'ensg' in expr.columns:
            expr.set_index('ensg', inplace=True)
        else:
            expr.set_index(expr.columns[0], inplace=True)

    met = [c for c in ('loc', 'name', 'tss_loc', 'strand') if c in expr.columns]
    if met:
        expr.drop(met, axis=1, inplace=True)
    inv = expr.columns[expr.dtypes == object].tolist()
    if inv:
        spacer.warning('\n')
        logger.warning('Invalid columns of datatype `object` (often text) '
                        'in expression data: {}\nThese columns will be '
                        'removed.'.format(inv))
        expr.drop(inv, axis=1, inplace=True)
    isna = expr.isna()
    if isna.any().any():
        spacer.error('\n')
        logger.error('Invalid expression data: data contains NaN values.')
        sys.exit(1)
    elif ctrl and (ctrl not in expr.columns.unique(0)):
        spacer.error('\n')
        logger.error('The control name of the samples `{}` was not found in '
                     'the passed expression data.'.format(ctrl))
        sys.exit(1)

    if expr.columns.nlevels > 1:
        exp_idx = [(name, dt) for name in expr.columns.unique(0) 
                    for dt in ['log2', 'z']]
        idx = expr.columns.values.tolist()
        misma = list(filter(lambda i: i not in exp_idx, idx))
        if any(misma):
            spacer.error('')
            msg = ('\tInvalid expresion data. When passing data with log2- and '
                   'z-data, the columns must be a MultiIndex in which level 0 '
                   'holds the names: [`name1`, ...] and level 1 the data types:'
                   ' [`log2`, `z`]. Expected column indices ({}):\n\t\t{}\n\t '
                   'Passed, unexpected column indices ({}):\n\t\t{}'
                    .format(len(exp_idx), exp_idx, len(misma), misma))
            logger.error(msg)
            sys.exit(1)
        else:
            return expr
    else:
        return util._add_log2_z(expr)

def _format_diff_genes(diff_genes_dir, genelists_mgtype='up', type_name=None): 
    """Take user gene list input input and format

    A single directory with deseq2 output files, a single dir. with up-genelist 
    files or 2 dirs. with up- and down- genelists are formatted here. A bool 
    DataFrame that holds the up- (and optionally down) differential genes is 
    returned.

    Args:
        diff_genes_dir: deseq2 directory, up-genelists dir. or list of up- and 
        down genelist dirs..
    genelists_mgtype: Which genelist type to handle. Internally used for 
        recursion 
    type_name: 'targets' or 'samples', depending on caller

    Returns:
        formatted _diff DataFrame
    """
    # check if path exists and contains TSV files
    def check_path(direc):
        if not os.path.exists(direc):
            spacer.info('')
            logger.error('Could not change directory to {}\nCheck the '
                        'path.'.format(os.path.abspath(direc)))
            sys.exit(1)
        files = glob.glob(direc + '/*.tsv')
        if not files:
            spacer.info('')
            logger.error('No *.tsv files found in {}\nCheck the path.'
                        .format(os.path.abspath(direc)))
            sys.exit(1)
            
    # check if up and down genelist directories are compatible
    def check_up_down_genelists():
        # check if 2 elements were passed, i.e. up+down genelist input
        if isinstance(diff_genes_dir, (list, tuple)) and len(diff_genes_dir) == 2:
            # check if paths are valid 
            check_path(diff_genes_dir[0])
            check_path(diff_genes_dir[1])
            # get the single TSV filenames
            up_dir = glob.glob(diff_genes_dir[0]+'/*.tsv')
            down_dir = glob.glob(diff_genes_dir[1]+'/*.tsv')
            # up and down must have the same number of elements
            if len(up_dir) != len(down_dir):
                msg = ('Number of up- and down genelist files differ. Found {} '
                       '*.tsv files in up directory\n{}\n{} *tsv files in down '
                       'directory:\n{}\n'.format(len(up_dir), 
                       os.path.abspath(diff_genes_dir[0]), len(down_dir), 
                       os.path.abspath(diff_genes_dir[1])))
                logger.error(msg)
                sys.exit(1)
            # to match up and down together safely, filenames must be the same
            f_up = [f[f.rfind(os.sep)+1:] for f in up_dir]
            f_down = [f[f.rfind(os.sep)+1:] for f in down_dir]
            is_single = lambda n: (n not in f_up) or (n not in f_down)
            singles = list(filter(is_single, set((*f_up, *f_down))))
            if singles:
                logger.error('Names of up- and down genelist files differ. '
                             'Names only found in one of the two '
                             'directories ({}):\n{}'
                             .format(len(singles), singles))
                sys.exit(1)
            # return the the up directory and that down mgs were passed
            return diff_genes_dir[0], True
        # a list of len=1 is treated as one element, both don't have down mgs
        elif isinstance(diff_genes_dir, (list, tuple)):
            check_path(diff_genes_dir[0])
            return diff_genes_dir[0], False
        else:
            check_path(diff_genes_dir)
            return diff_genes_dir, False

    # check whether deseq2 files or genelist files were passed
    def check_input_type():
        test_df = pd.read_csv(files[0], sep='\t')
            # check test dataframe for proper ensg index
        index_col = 'ensg'
        if 'ensg' not in test_df.columns:
            index_col = test_df.columns[0]
            if not str(test_df[index_col][0]).startswith('ENS'):
                spacer.error('')
                logger.error('The *.tsv files holding the gene keys do not '
                             'have  a column `ENS*` nor do they seem to have '
                             'an ensg index in the first column: {}, {}, ...'
                            .format(*test_df[index_col][:2].tolist()))
                sys.exit(1)
        # deseq2 output is identified based on the column names
        deseq2_cols = ['Unnamed: 0', 'baseMean', 'log2FoldChange', 'lfcSE', 
                       'stat', 'pvalue', 'padj']
        inp_type =  'deseq2' if test_df.columns.tolist() == deseq2_cols \
                            else 'genelist ({})'.format(genelists_mgtype)
        return inp_type, index_col

    # glue files together; return a list of pd.Series (columns)
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
                # genelists_mgtype is `down` in recursive function call
                s = pd.Series(True, index=d.index, name=(genelists_mgtype, name))
                diffs.append(s)
        return diffs

    direc, get_down_genelists = check_up_down_genelists()
    files = glob.glob(direc + '/*.tsv')
    inp_t, index_col = check_input_type()
    # if the samples are initiated from genelists, down mgs are required
    if type_name == 'samples' and inp_t != 'deseq2' and not get_down_genelists:
        spacer.error('')
        logger.error('When initiateing the samples diff. genes from genelist '
                     'input, both an up- and down directory with respective '
                     'genelists must be passed.')
        sys.exit(1)

    spacer.info('')
    st_st = [f[f.rfind(os.sep)+1:] for f in (files[0], files[-1])]
    f_type = inp_t if inp_t == 'deseq2' else genelists_mgtype
    logger.info('Formatting differential genes from {} files. {} *.tsv files '
                'in {}:\n{} ... {}\n'.format(f_type, len(files), direc, *st_st))

    diffs = merge_files()
    if inp_t == 'genelist (up)' and get_down_genelists:
        # for genelist data with down mgs, run function recursively with a 
        # single directory input, the down directory 
        diffs.extend(_format_diff_genes(diff_genes_dir[1], 'down'))
    elif inp_t == 'genelist (down)':
        # inside recursion: exit
        return diffs
    return pd.concat(diffs, axis=1, sort=True).fillna(False).sort_index(axis=1)