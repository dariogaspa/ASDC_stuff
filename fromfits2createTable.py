#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 12:30:24 2022

@author: dariogasparrini
"""
import numpy as np
from astropy.table import Table
import pandas
from collections import OrderedDict
#import astrotable_to_sql


def create_data_dictionaries(d):
    """
    Takes an astropy table "d" and extract out the data and meta data
    (units, description and format)
    :param d: astropy table
    :return dictdata: dictionary, keys are column names, with data stored as
                      values for each key (column)
    :return metadata: dictionary, keys are "columns, formats, units, description
                      values are the corresponding data from the astropy table
    """
    formats = dict(f='DOUBLE', b='BOOLEAN', i='INT', S='VARCHAR')
    dictdata = OrderedDict()
    col_descs = OrderedDict()
    col_formats = OrderedDict()
    col_units = OrderedDict()
 #   col_size = OrderedDict()
    # iterate around columns to scrape data
    for col in d.colnames:
        dictdata[col] = np.array(d[col].data)
        col_descs[col] = d[col].description
        col_units[col] = d[col].unit
        
        rawformat = d[col].descr[1]
        sizeformat= rawformat[2:]
        for f in formats:
            if f in rawformat:
                 col_formats[col] = formats[f]
                 if f == "S":
                    col_formats[col] = formats[f]+"("+sizeformat+")"
        if col not in col_formats:
            col_formats[col] = "STRING"
        #print col, col_formats,f
    # convert to metadata dictionary
    metadata = OrderedDict()
    metadata['columns'] = list(d.colnames)
    metadata['formats'] = list(col_formats.values())
    metadata['units'] = list(col_units.values())
    metadata['description'] = list(col_descs.values())
    # return dictionary of data and metadata
    return dictdata, metadata

# =============================================================================
# Define variables
# =============================================================================
WORKSPACE = '/Users/dariogasparrini/Downloads/'
# file fits
DATAPATH = WORKSPACE + 'gll_psc_v30_exp.fit'
# Name the sql table the catalogue should create
CAT_NAME = 'mmc_fermi4fgldr3'
# -----------------------------------------------------------------------------
# Load data with Astropy table
print("\n Loading table with astropy (slow)...")
data = Table.read(DATAPATH)
# -------------------------------------------------------------------------
# Convert to pandas dataframe via ordered dictionary
print("\n Loading data into pandas data frame...")
ddata, dmeta = create_data_dictionaries(data)
pandas_data = pandas.DataFrame(ddata)
pandas_data = pandas_data.replace([np.inf, -np.inf], np.nan)
# -------------------------------------------------------------------------
# Create meta data pandas dataframe
print("\n Loading meta data into pandas data frame...")
pandas_meta = pandas.DataFrame(dmeta)
pandas_meta = pandas_meta.replace([np.inf, -np.inf], np.nan)
print "CREATE TABLE "+ CAT_NAME + "(\n"
with pandas.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print pandas_meta.iloc[:,[0,1]]
print ")\n"
# -------------------------------------------------------------------------
