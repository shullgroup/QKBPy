# -*- coding: utf-8 -*-
"""
Created on Wed 13 05 2026 20:35:00
@author: brodericklewis
"""
import csv
from cycler import cycler
#plt.rcParams['axes.prop_cycle']
default_cycler = cycler(color=[
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
    except ValueError:
        return False
    
def first_line(path, **kwargs):
    '''
    Find first line of data file to skip headers.
    
    Parameters
    ----------
    path : Path
        Path to data file using pathlib Path object.

    delim : str, default '\t'
        Delimiter used for reading the file.
    
    Returns
    -------
    start_row : int
        First row in file with numeric data.

    '''

    delim = kwargs.get('delim', '\t')
    
    with open(path, 'r') as f:
            #check which row starts the data
            reader = csv.reader(f, delimiter=delim)
            for i, row in enumerate(reader):
                if any(is_numeric(cell) for cell in row if cell.strip()):
                    start_row = i
                    break

    return start_row

def first_line(path, **kwargs):
    '''
    Find first line of data file to skip headers.
    
    Parameters
    ----------
    path : Path
        Path to data file using pathlib Path object.
    sep : str, default ','
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
    sep = kwargs.get('sep', ',')
    target_cols = kwargs.get('target_cols', None)
    encoding = kwargs.get('encoding', 'utf-8')
    
    # Try opening with UTF-8, fallback to latin-1 for special characters like 
    try:
        with open(path, 'r', encoding=encoding) as f:
            f.readline()
    except UnicodeDecodeError:
        encoding = 'latin-1'

    with open(path, 'r', encoding=encoding) as f:
        reader = csv.reader(f, delimiter=sep)
        for i, row in enumerate(reader):
            # Strip whitespace from cells
            cleaned_row = [cell.strip() for cell in row]
            
            # Target specific columns
            if target_cols is not None:
                if len(cleaned_row) > max(target_cols):
                    # Check if all columns in our target list are numeric
                    if all(is_numeric(cleaned_row[c]) for c in target_cols):
                        return i

    return 0
