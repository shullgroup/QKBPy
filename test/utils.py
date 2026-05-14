# -*- coding: utf-8 -*-
"""
Created on Wed 13 05 2026 20:35:00
@author: brodericklewis
"""
import csv
import pandas as pd
from cycler import cycler

# Define the cycler here as a constant
DEFAULT_CYCLER = cycler(color=[
    '#0093F5', '#F08E2C', '#000000', '#424EBD', '#B04D25', '#75CA85', '#C892D6'
]*3, linestyle=['-']*7 + ['--']*7 + [':']*7)

def is_numeric(val):
    try:
        float(val)
        return True
    except (ValueError, TypeError):
        return False

def find_data_start(path, sep=',', target_cols=None, encoding='utf-8'):
    """
    Finds the first row index containing numeric data.
    Standardized to handle encoding fallbacks.
    """
    try:
        with open(path, 'r', encoding=encoding) as f:
            f.readline()
    except UnicodeDecodeError:
        encoding = 'latin-1'

    with open(path, 'r', encoding=encoding) as f:
        reader = csv.reader(f, delimiter=sep)
        for i, row in enumerate(reader):
            cleaned_row = [cell.strip() for cell in row if cell.strip()]
            if not cleaned_row:
                continue
            
            # If target_cols is provided, check only those
            if target_cols:
                if len(cleaned_row) > max(target_cols):
                    if all(is_numeric(cleaned_row[c]) for c in target_cols):
                        return i
            # Otherwise, check if the majority of the row is numeric
            else:
                if sum(is_numeric(c) for c in cleaned_row) >= 2:
                    return i
    return 0

def apply_lib_style(ax):
    """Utility to apply the library's cycler to an axis."""
    ax.set_prop_cycle(DEFAULT_CYCLER)
