# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 13:52:21 2024

@author: Collin
"""

#### reads in all files within a directory, containing .txt files with user inputs for defining variables ####
#### reads in a dataset and calculates a product specified by user input ####

def plot_reg_ann_data(inputs_path):
    
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from funcs.make_vars import make_vars
    from funcs.make_poly import make_poly
    from funcs.make_poly_mask import make_poly_mask
    from funcs.limit_data_bbox import limit_data_bbox
    from funcs.calc_Z import calc_Z
    from funcs.get_shp_coords import get_shp_coords
    
    # Read user inputs from the .txt file
    with open(inputs_path, 'r') as file:
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
    else:
        d = d

    # calls make_poly function that creates polygon and arrays of corresponding coords 
    xb,yb,basin,poly,xint,yint,polyint,polyc = make_poly(reg_choice,use_polygon_coast)

    # if only looking at specified region, create mask that masks pts outside region
    if perform_mask_calculation.lower() == 'y':
        x_poly, y_poly, mask = make_poly_mask(d, polyc, x_grid, y_grid)
    else:
        mask = np.zeros_like(d_ann, dtype=bool)
    
    # define Bamber DEM coordinates
    xmin_bamber=-560.*5e3
    xmax_bamber= 560.*5e3
    ymax_bamber= 560.*5e3
    ymin_bamber=-560.*5e3

    # define x,y grid spacing
    dx=grid_res*1000
    dy=grid_res*1000

    x_bamber = np.arange(xmin_bamber,xmax_bamber+dx,dx)
    y_bamber = np.arange(ymax_bamber,ymin_bamber-dy,-dy)
    
    # limit plotting coordinates by bounding box according to region 
    x_min, x_max, y_min, y_max = limit_data_bbox(x_bamber, y_bamber, polyc, 1.5e6, 1.5e6)

    # use calc_Z to create masked data array using region polygon    
    d_mask, ref, Z_mon = calc_Z(d,mask,z0,zf)
    
    # create annual data array masked by region polygon   
    d_ann_mask = np.ma.array(d_ann, mask=mask[:len(d_ann)])
    
    # calculate Z-scores for region data relative to z0 and zf, defined by input file 
    Z = (d_ann_mask - np.mean(d_ann_mask[z0:zf,:,:], axis=0)) / np.std(d_ann_mask[z0:zf,:,:], axis=0)

    start_yr = 1979
    if which_set == 1:
        start_yr = 1980
    elif which_set == 3:
        start_yr = 1980
        
    variable = 'z'
    z_option = 'dif'
      
    #variable = input("Enter 'z' to plot z scores or '{}' to plot annual {}, Moving Average, Moving Standard Deviation, or Replace {}th+ w/ normal: ".format(dtype,dtype,perc_repl_thresh))

    if variable.lower() == 'z':
        
        #z_option = input("Enter 'original' to plot Z or 'dif' to plot Z_dif: ")
        
        if z_option.lower() == 'original':
            dplt = Z
            dtag = 'Z'
            levels = np.arange(-4, 4.1, 0.5)  # Levels for z scores
            label = 'Z Scores'
            cmap = 'RdBu'
        elif z_option.lower() == 'dif':
            from funcs.calc_repl_normal import calc_repl_normal
            d90, d_norm_repl, d_ann_norm_repl = calc_repl_normal(d_mask, ref, perc_repl_thresh, 100)
            Z_REPL = (d_ann_norm_repl - np.mean(d_ann_mask[z0:zf,:,:], axis=0)) / np.std(d_ann_mask[z0:zf,:,:], axis=0)
            dplt = Z - Z_REPL
            dtag = 'Z_DIF'
            levels = np.arange(-4, 4.1, 0.5)
            label = 'Difference in Z score (SMB vs 90th perc replace normal vals)'
            cmap = 'RdBu'
    elif variable.lower() == dtype:
        data_variable = input("Enter 'ann' to plot ANN SMB, 'ma' to plot Moving Average, 'sd' to plot Moving Standard Deviation, or 'norm' for Replace 90th+ w/ normal': ")
        if data_variable.lower() == 'ann':
            dplt = d_ann_mask
            dtag = 'ann'
            #levels = np.linspace(0, 250, 10)
            levels = np.linspace(0, 300, 11)
            #levels = np.arange(np.nanpercentile(data, 68.), np.nanpercentile(data, 99.7), (np.nanpercentile(data, 99.7) - np.nanpercentile(data, 68)) / 10)  # Levels for ANN SMB
            label = 'Annual {}'.format(dtype)
            cmap = 'viridis'
        elif data_variable.lower() == 'ma':
            from funcs.calc_moving_avg import calc_moving_avg
            ma, msd = calc_moving_avg(d_ann_mask, zf-z0)
            # Calculate moving average
            dplt = ma
            dtag = '{}yr_ma'.format(str(zf-z0))
            levels = np.arange(np.nanpercentile(dplt, 68), np.nanpercentile(dplt, 96), (np.nanpercentile(dplt, 96) - np.nanpercentile(dplt, 68)) / 10)  # Levels for moving average
            #levels = np.arange(-200, 1000, 100)
            label = 'Moving Average'
            cmap = 'viridis'
        elif data_variable.lower() == 'sd':
            from funcs.calc_moving_avg import calc_moving_avg
            ma, msd = calc_moving_avg(d_ann_mask, zf-z0)
            # Calculate moving standard deviation
            dplt = msd
            dtag = '{}yr_msd'.format(str(zf-z0))
            levels = np.arange(np.percentile(dplt, 0.02), np.percentile(dplt, 99.), (np.percentile(dplt, 99.) - np.percentile(dplt, 0.1)) / 10)  # Levels for moving standard deviation
            #levels = np.arange(-200, 1000, 10)  # Levels for moving average
            label = 'Moving Standard Deviation'
            cmap = 'viridis'
        elif data_variable.lower() == 'norm':
            from funcs.calc_repl_normal import calc_repl_normal
            d90, d_norm_repl, d_ann_norm_repl = calc_repl_normal(d_ann_mask, ref, perc_repl_thresh,100)
            dplt = d_ann_norm_repl 
            dtag = 'norm_repl'
            levels = np.arange(np.nanpercentile(dplt, 68), np.nanpercentile(dplt, 99), (np.nanpercentile(dplt, 99) - np.nanpercentile(dplt, 68)) / 10)  # Levels for moving average
            #levels = np.arange(0, 1000, 50)
            label = 'Replace 90th+ w/ normal'
            cmap = 'Blues'
        
        else:
            print("Invalid input. Please try again.")
            return
    else:
        print("Invalid input. Please try again.")
        return
    
    polygons = get_shp_coords('ANT')
    
    for yr in range(0, dplt.shape[0]):
        plt.figure(figsize=(10,8))
    
        if poly is not None:
            plt.plot(x_poly, y_poly, 'k')
            plt.plot(xint, yint, color='m')
    
        for polygon in polygons:
            xi,yi = polygon
            plt.plot(xi, yi,color='k',lw=1.5, zorder=4)
        plt.plot(xint, yint, 'm', lw = 2.5)
        
        res = plt.contourf(x_grid, y_grid, dplt[yr, :, :], levels=levels, cmap=cmap, extend='both',zorder=1)

    # Add contour lines
        contours = plt.contour(x_grid, y_grid, dplt[yr, :, :], levels=levels, colors='k', linewidths=0.35)
        plt.clabel(contours, inline=True, fontsize=7)

        plt.xlim(x_min, x_max)
        plt.ylim(y_min, y_max)
        
        #plt.xlim(-.9e6, 2.35e6)
        #plt.ylim(.75e6, 2.25e6)
        
        #x_min, x_max = np.min(x_grid), np.max(x_grid)
        #y_min, y_max = np.min(y_grid), np.max(y_grid)
        
        #plt.xlim(-0.9e6, x_max)
        #plt.ylim(y_min+0.5e6, y_max-0.5e6)
        
        plt.grid()
        plt.colorbar(res, orientation='vertical', label=label, extend='both', fraction=0.046, pad=0.04)
        
        plt.gca().set_xticks([])  # Remove x-axis tick labels
        plt.gca().set_yticks([])  # Remove y-axis tick labels

        plt.title('{}'.format(start_yr + yr))
 
        plt.savefig('D:/Research/PRODUCTS/reg_ann/cesm2/{}_{}_ann_{}_{}_{}_{}km.jpg'.format(dset, dtype, dtag, start_yr + yr, reg_choice, grid_res))
        
        plt.close()

import os
dir_path = 'D:/Research/SCRIPTS/crich/user_input/plot_reg_ann/cesm2_smb/27km/'
for f_path in os.listdir(dir_path): 
    f = os.path.join(dir_path, f_path)
    plot_reg_ann_data(f)