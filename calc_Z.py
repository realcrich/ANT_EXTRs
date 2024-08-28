# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 18:44:08 2024

@author: Collin
"""

# calculate Z-score, data array masked by region polygon, and ref (mean value over reference period)

# inputs: 
    # data - yep, data
    # mask - the region mask calculated from polygon for region of interest
    # ref_start - index of first year to calculate Z-scores relative to 
    # ref_end - "       ... last year ...                           "

import numpy as np
def calc_Z(data, mask, ref_start, ref_end):

    # mask the data array using region polygon
    data_mask = np.ma.array(data, mask=mask)

    # Calculate the reference mean using the specified years
    ref = np.ma.mean(data_mask[ref_start:ref_end], axis=0)

    # Calculate the Z score 
    z = (data_mask - ref) / np.ma.std(data_mask[ref_start:ref_end, :, :], axis=0)

    return data_mask, ref, z 