# -*- coding: utf-8 -*-
# updated version of this file is maintained at
# https://github.com/shullgroup/QKBPy/blob/main/test/utils.py

import numpy as np
import csv
from cycler import cycler

# Shared across the library
DEFAULT_CYCLER = cycler(color=[
    '#0093F5', '#F08E2C', '#000000', '#424EBD', '#B04D25', '#75CA85', '#C892D6'
]*3, linestyle=['-']*7 + ['--']*7 + [':']*7)

def is_numeric(cell):
    '''
    Check if a cell has numeric data. Used to start reading data files.
    
    Parameters
    ----------
    cell : str
        String which may be numeric (cell in a given row)
    
    Returns
    -------
    bool
        If cell is numeric (True) or not (False)
    '''
    try:
        float(cell)
        return True
    except (ValueError, TypeError):
        return False

def first_line(path, **kwargs):
    '''
    Find first line of data file to skip headers.
    
    Parameters
    ----------
    path : Path
        Path to data file using pathlib Path object.
    sep : str, default '\t'
        Delimiter used in file.
    target_cols : list, default None
        Array/list of specific columns to read.
    encoding : str, default 'utf-8'
        Encoding style option if this is ever a problem.
    
    Returns
    -------
    start_row : int
        First row in file with numeric data.
    '''
    sep = kwargs.get('sep', '\t')
    target_cols = kwargs.get('target_cols', None)
    encoding = kwargs.get('encoding', 'utf-8')
    
    # Try opening with UTF-8, fallback to latin-1
    try:
        with open(path, 'r', encoding=encoding) as f:
            f.readline()
    except UnicodeDecodeError:
        encoding = 'latin-1'

    with open(path, 'r', encoding=encoding) as f:
        reader = csv.reader(f, delimiter=sep)
        for i, row in enumerate(reader):
            cleaned_row = [cell.strip() for cell in row]
            if not cleaned_row:
                continue
            
            if target_cols is not None:
                if len(cleaned_row) > max(target_cols):
                    if all(is_numeric(cleaned_row[c]) for c in target_cols):
                        return i
            else:
                if any(is_numeric(cell) for cell in cleaned_row):
                    return i
    return 0

def remove_step_lines(df):
    '''
    Remove extra lines in files from multiple steps in 
    TA Instruments experiments
    
    Parameters
    ----------
    df : pd.DataFrame
        Dataframe from experimental file
    
    Returns
    -------
    df : pd.DataFrame
        Cleaned dataframe with step lines removed.
    '''

    to_drop = []
    for rows in np.arange(0,len(df)):
        
        if df.iloc[rows, 0] == '[step]':
            to_drop.extend([rows,rows+1,rows+2,rows+3])
    
    df = df.drop(to_drop).reset_index(drop=True)

    return df