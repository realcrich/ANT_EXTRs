# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#### reads in all files within a directory, containing .txt files with user inputs for defining variables ####

import os
import numpy as np
import matplotlib.pyplot as plt
from funcs.make_vars import make_vars
from funcs.make_poly import make_poly
from funcs.make_poly_mask import make_poly_mask
from funcs.calc_repl_normal import calc_repl_normal
from funcs.calc_Z import calc_Z
from funcs.calc_moving_avg import calc_moving_avg

def plot_MA_time_series_percs(file_path):
    
    # Read user inputs from the .txt file
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Extract inputs from lines
    base_dir = lines[0].strip()
    in_file = os.path.join(base_dir, lines[1].strip())
    rem_ref = lines[2].strip()
    which_set = int(lines[3].strip())
    dtype = lines[4].strip()
    reg_choice = lines[5].strip()
    use_polygon_coast = lines[6].strip()
    perform_mask_calculation = lines[7].strip()
    grid_res = int(lines[8].strip())
    perc_repl_thresh = int(lines[9].strip())
    z0 = int(lines[10].strip())
    zf = int(lines[11].strip())
    
    # create tag for dataset for plotting functions
    if which_set == 0: 
        dset = 'RACMO2.3p2'
    elif which_set == 1: 
        dset = 'MERRA2'
    elif which_set == 2: 
        dset = 'CESM2'
    elif which_set == 3: 
        dset = 'UKESM'
    else: 
        print('Invalid entry for dataset')
    
    # Call make_vars function with in_file
    d, d_ann, time, lon, lat, x_grid, y_grid = make_vars(in_file,dtype,grid_res)

    # For ANT reference period is historical - 2008, due to anomalous accumulation in QML
    # For RACMO2.3p2 and CESM2 historical start is 1979 and for MERRA2 and UKESM, 1980
    # aka are we looking at anomalies or no?
    if np.logical_and(rem_ref == 'y', which_set == '0'):
        d = d - np.mean(d[0:348,:,:], axis=0)
    elif np.logical_and(rem_ref == 'y', which_set == '1'):
        d = d - np.mean(d[0:336,:,:], axis=0)
    elif np.logical_and(rem_ref == 'y', which_set == '2'):
        d = d - np.mean(d[0:348,:,:], axis=0)
    elif np.logical_and(rem_ref == 'y', which_set == '3'):
        d = d - np.mean(d[0:336,:,:], axis=0)
    elif rem_ref == 'n':
        d = d

    # calls make_poly function that creates polygon and arrays of corresponding coords
    xb,yb,basin,poly,xint,yint,polyint,polyc = make_poly(reg_choice,use_polygon_coast)

    # if only looking at specified region, create mask that masks pts outside region
    if perform_mask_calculation.lower() == 'y':
        x_poly, y_poly, mask = make_poly_mask(d, polyc, x_grid, y_grid)
    else:
        mask = np.zeros_like(d, dtype=bool)
        _, _, mask = None, None, None

    # use calc_Z to create masked data array using region polygon
    d_mask, ref, Z = calc_Z(d,mask,z0,zf)
    
    d_ann_mask,_,_ = calc_Z(d_ann, mask[:len(d_ann)], z0, zf)

    # create masks for all cells within percentile ranges, as well as subsituting normal replace
    d70, d70_norm_repl, d70_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, 70, 89)
    d90, d90_norm_repl, d90_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, 90, 98)
    d99, d99_norm_repl, d99_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, 99, 100)

    # masks region data to only within percentile ranges
    d70_mask = np.ma.array(d_mask, mask = d70)
    d90_mask = np.ma.array(d_mask, mask = d90)
    d99_mask = np.ma.array(d_mask, mask = d99)
    
    nmons = len(d)
    nyrs = nmons//12
    int_yrs = (nmons//12)*12
    
    d70_mask_yr = np.ma.sum(d70_mask[:int_yrs].reshape(d70.shape[0]//12, 12, d70_mask.shape[1], d70_mask.shape[2]),axis=1)
    d90_mask_yr = np.ma.sum(d90_mask[:int_yrs].reshape(d90.shape[0]//12, 12, d90_mask.shape[1], d90_mask.shape[2]),axis=1)
    d99_mask_yr = np.ma.sum(d99_mask[:int_yrs].reshape(d99.shape[0]//12, 12, d99_mask.shape[1], d99_mask.shape[2]),axis=1)
    
    # calculate portion of annual accumlation that comes from percentile ranges
    yr_cont70 = np.ma.sum(d70_mask_yr, axis=(1,2)) / np.ma.sum(d_ann_mask, axis=(1,2))
    yr_cont90 = np.ma.sum(d90_mask_yr, axis=(1,2)) / np.ma.sum(d_ann_mask, axis=(1,2))
    yr_cont99 = np.ma.sum(d99_mask_yr, axis=(1,2)) / np.ma.sum(d_ann_mask, axis=(1,2))

    # calculate portion of annual accumlation that comes from percentile ranges
    mon_cont70 = np.ma.sum(d70_mask, axis=(1,2)) / np.ma.sum(d_mask, axis=(1,2))
    mon_cont90 = np.ma.sum(d90_mask, axis=(1,2)) / np.ma.sum(d_mask, axis=(1,2))
    mon_cont99 = np.ma.sum(d99_mask, axis=(1,2)) / np.ma.sum(d_mask, axis=(1,2))

    # calculate moving average on percentage arrays
    ma70, msd70 = calc_moving_avg(yr_cont70, 5)
    ma90, msd90 = calc_moving_avg(yr_cont90, 5)
    ma99, msd99 = calc_moving_avg(yr_cont99, 5)

    # set start year based on dset
    start_yr = 1979
    if which_set == 1:
        start_yr = 1980
    elif which_set == 3:
        start_yr = 1980
    nmons = len(d)   
    end_yr = start_yr + (nmons/12)

    plt.figure(figsize=(8,5))
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), yr_cont70, marker='o', ms = 4, label='70-89th') 
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), yr_cont90, marker='+', ms = 4, label='90-98th')
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), yr_cont99, marker='d', ms = 4, label='99th+')
    plt.legend()
    plt.ylabel('Contribution to Annual {} from Percentile Values'.format(dtype))
    plt.savefig('D:/Research/PRODUCTS/percentiles/no_smooth/{}_{}_EPE_classes_{}_{}km.png'.format(dset,dtype,reg_choice,grid_res))
    plt.close()

# define directory path where input files are located on device
# looping thru input files for each region , run plot_distribution 
# save .png files to specified location on device    
#dir_path = 'D:/Research/SCRIPTS/crich/user_input/Z/racmo_smb/27km/'
dir_path = 'D:/Research/SCRIPTS/crich/user_input/racmo_updated_timeseries/'
for f_path in os.listdir(dir_path): 
    f = os.path.join(dir_path, f_path)
    plot_MA_time_series_percs(f)