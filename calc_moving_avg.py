# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 11:09:14 2024

@author: Collin
"""

# calculates a moving average of input data array using a user specified window
# size, if the window size is not an integer multiple of the size of the data array 
# the window size is shortened to half and stepped back up to the whole window size 
# and at the beginning and end of the data array 

import numpy as np

def calc_moving_avg(data, window_size):
    
    # cut window size in half to use for stepping at beginning/end of array 
    window_half = window_size // 2
    
    # create zero value arrays for the mean and standard deviations to be calculated
    moving_avg = np.zeros_like(data)
    moving_sd = np.zeros_like(data)
    
    # start index
    yr = 0 
    start = 0
    end = window_half
    nyrs = len(data)

    # calculate averages over indices with window size adjusting according to len(data)
    while yr <= nyrs - 1:

        if yr == 0:
            window = data[start : end]
            end+=1

        elif yr <= window_half:
            window = data[start : end]
            end+=1
        elif yr >= nyrs - window_half - 1: 
            window = data[start:end]
            start+=1

        elif np.logical_or((yr <= nyrs - 1), (yr >= window_half)):
            window = data[yr - window_size // 2 : yr + window_size // 2]
            start+=1
            end+=1

        # for each year calculate the mean value over the window and increment up 
        moving_avg[yr] = np.mean(window, axis=0)
        moving_sd[yr] = np.std(window, axis=0)
        yr+=1


    return moving_avg, moving_sd