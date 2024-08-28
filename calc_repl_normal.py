# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 11:29:47 2024

@author: Collin
"""

# take a copy of input data array and replace values above a user specified threshold 
# with the mean value for that grid cell over the reference period 

# inputs: 
    # data - explains itself
    # ref - the grid of mean values calculated over the reference period
    # thresh_bot - lower threshold for percentile range
    # thresh_top - upper threshold for percentile range

import numpy as np

def calc_repl_normal(data, ref, mask, thresh_bot, thresh_top):
    
    # arrays for the values associated with lower and upper percentile thresholds
    perc_bot = np.percentile(data, thresh_bot, axis=0)
    perc_top = np.percentile(data, thresh_top, axis=0)

    # create mask array where all values are set to boolean False
    perc_mask = np.ones(data.shape,dtype=bool)

    # set grid cells where the value is both below the upper threshold 
    # AND above the lower threshold to False (we are interested in looking at these values) 
    perc_mask = np.where(np.logical_and(np.logical_and(data>=perc_bot,data<=perc_top), mask==False), False, True)

    # Use the mask to replace values in data_norm_repl with ref
    data_norm_repl = np.where(perc_mask, ref, data)
    
    num_months = data_norm_repl.shape[0]
    num_years = num_months // 12

    # Reshape data_norm_repl to have three dimensions (years, months_per_year, lon, lat)
    reshaped_data = data_norm_repl[:num_years * 12].reshape(num_years, 12, *data_norm_repl.shape[1:])

    # Calculate the sum along the months_per_year axis
    data_ann_norm_repl = np.sum(reshaped_data, axis=1)

    return perc_mask, data_norm_repl, data_ann_norm_repl