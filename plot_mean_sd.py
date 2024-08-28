# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:15:07 2024

@author: Collin
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from funcs.make_vars import make_vars
from funcs.make_poly import make_poly
from funcs.make_poly_mask import make_poly_mask
from funcs.calc_repl_normal import calc_repl_normal
from funcs.calc_Z import calc_Z
#from funcs.calc_moving_avg import calc_moving_avg

def calc_moving_avg(data, window_size):

    # Create zero value arrays for the mean and standard deviations to be calculated
    moving_avg = np.zeros(data.shape[0])
    moving_sd = np.zeros(data.shape[0])

    nyrs, nlon, nlat = data.shape

    # Calculate moving average and standard deviation along the first axis (nyrs)
    for yr in range(nyrs):
        start = max(0, yr - window_size // 2)
        end = min(nyrs, yr + window_size // 2 + 1)
        window = data[start:end]

        # Calculate mean and standard deviation over the window along the basin's dimensions
        basin_avg = np.mean(window, axis=(1, 2))
        basin_sd = np.std(window, axis=(1, 2))

        # Store the basin-wide average and standard deviation for the current year
        moving_avg[yr] = np.mean(basin_avg)
        moving_sd[yr] = np.mean(basin_sd)

    return moving_avg, moving_sd

def plot_mean_sd(file_path):
    
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
        
    if np.logical_and(rem_ref == 'y', which_set == '0'):
        d_ann = d_ann - np.mean(d_ann[0:29,:,:], axis=0)
    elif np.logical_and(rem_ref == 'y', which_set == '1'):
        d_ann = d_ann - np.mean(d_ann[0:28,:,:], axis=0)
    elif np.logical_and(rem_ref == 'y', which_set == '2'):
        d_ann = d_ann - np.mean(d_ann[0:29,:,:], axis=0)
    elif np.logical_and(rem_ref == 'y', which_set == '3'):
        d_ann = d_ann - np.mean(d_ann[0:28,:,:], axis=0)
    elif rem_ref == 'n':
        d_ann = d_ann

    from shapely.ops import cascaded_union
    # calls make_poly function that creates polygon and arrays of corresponding coords
    xb,yb,basin,poly,xint,yint,polyint,polyc = make_poly(reg_choice,use_polygon_coast)
    #_,_,_,_,_,_,_,polyapb = make_poly('apb',use_polygon_coast)
    #_,_,_,_,_,_,_,polyka = make_poly('ka',use_polygon_coast)
    
    #polyc = cascaded_union([polyc, polyapb, polyka])

    # if only looking at specified region, create mask that masks pts outside region
    if perform_mask_calculation.lower() == 'y':
        x_poly, y_poly, mask = make_poly_mask(d, polyc, x_grid, y_grid)
    else:
        mask = np.zeros_like(d, dtype=bool)
        _, _, mask = None, None, None

    # use calc_Z to create masked data array using region polygon
    d_mask, ref, Z = calc_Z(d,mask,z0,zf)
    
    d_ann_mask,_,_ = calc_Z(d_ann, mask[:len(d_ann)], z0, zf)
    
    nmons = len(d)
    nyrs = nmons//12
    int_yrs = (nmons//12)*12

    # set start year based on dset
    start_yr = 1979
    if which_set == 1:
        start_yr = 1980
    elif which_set == 3:
        start_yr = 1980
    nmons = len(d)   
    end_yr = start_yr + (nmons/12)
    
    # calculate moving average on percentage arrays
    dam_ma, dam_sd = calc_moving_avg(d_ann_mask, 10)

    plt.figure(figsize=(11,5))
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), dam_ma, color='b', marker='o', ms = 5, lw=1, label='ann mean SMB') 
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), dam_sd, color='orange',marker='+', ms = 5, lw=1, label='ann std SMB')
    plt.legend()
    plt.ylabel('SMB (mm w.e.)')
    plt.title('{}'.format(reg_choice))
    plt.savefig('D:/Research/PRODUCTS/MEAN_SD/{}_{}_{}_{}km.png'.format(dset,dtype,reg_choice,grid_res))
    plt.close()
    
# define directory path where input files are located on device
# looping thru input files for each region , run plot_distribution 
# save .png files to specified location on device    
dir_path = 'D:/Research/SCRIPTS/crich/user_input/Z/merra2_precip/27km/'
for f_path in os.listdir(dir_path): 
    f = os.path.join(dir_path, f_path)
    plot_mean_sd(f)