# -*- coding: utf-8 -*-
"""
Created on Wed Feb  7 10:33:39 2024

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
#from shapely.geometry import MultiPolygon
from shapely.ops import cascaded_union

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
    
    # set start year based on dset
    start_yr = 1979
    if which_set == 1:
        start_yr = 1980
    elif which_set == 3:
        start_yr = 1980
    nmons = len(d)   
    end_yr = start_yr + (nmons/12)
    
    nlon = len(lon)
    nlat = len(lat)

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
    #xb,yb,basin,poly,xint,yint,polyint,polyc = make_poly(reg_choice,use_polygon_coast)
    xb,yb,basin,poly,xint,yint,polyint,polyc = make_poly('aap',use_polygon_coast)
    #_,_,_,_,_,_,_,polyapb = make_poly('apb',use_polygon_coast)
    #_,_,_,_,_,_,_,polyka = make_poly('ka',use_polygon_coast)
    
    #polyc = cascaded_union([polyc, polyapb, polyka])
    
    # if only looking at specified region, create mask that masks pts outside region
    if perform_mask_calculation.lower() == 'y':
        x_poly, y_poly, mask = make_poly_mask(d, polyc, x_grid, y_grid)
    else:
        mask = np.zeros_like(d, dtype=bool)
        _, _, mask = None, None, None
        
    # limit plotting coordinates by bounding box according to region 
    #x_min, x_max, y_min, y_max = limit_data_bbox(x_bamber, y_bamber, polyc, 1.5e6, 1.5e6)
    zf = 30
    # use calc_Z to create masked data array using region polygon
    d_mask, ref, Z = calc_Z(d,mask,z0,zf)
    
    #d_ann_mask,_,_ = calc_Z(d_ann, mask[:len(d_ann)], z0, zf)
    d_ann_mask,ref,Z = calc_Z(d_ann, mask[:len(d_ann)], z0, zf)
    
    n_reg_pts = np.sum(~mask[0], axis=(0,1))
    #n_all_pts = np.sum(np.ones_like(d)[0],axis=(0,1))
    #npts = n_all_pts - n_reg_pts
    
    # create masks for all cells within percentile ranges, as well as subsituting normal replace
    d90, d90_norm_repl, d90_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, 90, 100)
    d98, d98_norm_repl, d98_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, 95, 100)
    d99, d99_norm_repl, d99_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, 99, 100)

    # masks region data to only within percentile ranges
    d90_mask = np.ma.array(d_mask, mask = d90)
    d98_mask = np.ma.array(d_mask, mask = d98)
    d99_mask = np.ma.array(d_mask, mask = d99)
    
    # Here we initiate an array of zeros for each class of EPE, wherever 
    # the threshold is met a value of 1 is given
    nEPE = np.zeros((3,nmons,nlon,nlat),dtype=int)
    nEPE[0] = np.where(d90 == False, 1, 0)
    nEPE[1] = np.where(d98 == False, 1, 0)
    nEPE[2] = np.where(d99 == False, 1, 0)
    #nEPE[0, ~d90.mask] = 1
    #nEPE[1, ~d98.mask] = 1
    #nEPE[2, ~d99.mask] = 1
    
    nyrs = nmons//12
    int_yrs = (nmons//12)*12
    
    nEPE_yr = np.sum(nEPE[:,:int_yrs,:,:].reshape(3, nyrs, 12, nlon, nlat), axis=2) / 12
    
    # sum over the lon, lat indices to find the number of cells per month that experience
    # an EPE over a threshold and divide by the total area * mons per yr to calculate ratio
    n90_yr = np.sum(nEPE_yr[0],axis=(1,2)) / (n_reg_pts)
    n98_yr = np.sum(nEPE_yr[1],axis=(1,2)) / (n_reg_pts)
    n99_yr = np.sum(nEPE_yr[2],axis=(1,2)) / (n_reg_pts)
    
    d90_ann_mask = np.sum(d90_mask[:int_yrs].reshape(nyrs, 12, nlon, nlat), axis=1)
    d95_ann_mask = np.sum(d98_mask[:int_yrs].reshape(nyrs, 12, nlon, nlat), axis=1)
    d99_ann_mask = np.sum(d99_mask[:int_yrs].reshape(nyrs, 12, nlon, nlat), axis=1)

    d90_mcont = np.sum(d90_ann_mask,axis=(1,2)) / np.sum(d_ann_mask,axis=(1,2))
    d95_mcont = np.sum(d95_ann_mask,axis=(1,2)) / np.sum(d_ann_mask,axis=(1,2))
    d99_mcont = np.sum(d99_ann_mask,axis=(1,2)) / np.sum(d_ann_mask,axis=(1,2))
    
    plt.figure(figsize=(14,6))
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), n90_yr, 'b', marker='o', ms = 4, lw = 1, ls = '--', label='90th+ freq') 
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), n98_yr, 'r',marker='+', ms = 4, lw = 1, ls = '--',label='95th+ freq')
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), n99_yr, 'orange', marker='d', ms = 4, lw = 1, ls = '--',label='99th+ freq')

    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), d90_mcont, 'b', lw = 6, label='90th+ m cont', alpha=0.3) 
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), d95_mcont, 'r', lw = 6, label='95th+ m cont', alpha=0.3)
    plt.plot(np.linspace(start_yr, end_yr+1, nyrs), d99_mcont, 'orange', lw = 6, label='99th+ m cont', alpha=0.3)

    plt.legend(fontsize='8')
    plt.title('{}'.format('DML'))
    plt.ylabel('ratio of (nEPEs/npts) & (mEPEs/ann smb)')
    #plt.savefig('D:/Research/PRODUCTS/nEPEs/{}_{}_EPE_90_98_99_{}_{}km_QML.png'.format(dset,dtype,reg_choice,grid_res))
    plt.savefig('D:/Research/PRODUCTS/nEPEs/{}_{}_EPE_90_95_99_{}km_AAP.png'.format(dset,dtype,grid_res))
    plt.close()
    
# define directory path where input files are located on device
# looping thru input files for each region , run plot_distribution 
# save .png files to specified location on device    
#dir_path = 'D:/Research/SCRIPTS/crich/user_input/Z/merra2_precip/27km/'
#dir_path = 'D:/Research/SCRIPTS/crich/user_input/racmo_updated_timeseries/'
dir_path = 'D:/Research/SCRIPTS/crich/user_input/racmo2.4/'
for f_path in os.listdir(dir_path): 
    f = os.path.join(dir_path, f_path)
    plot_MA_time_series_percs(f)       
    