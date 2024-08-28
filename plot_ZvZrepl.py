# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 12:05:48 2024

@author: Collin
"""
#### reads in all files within a directory, containing .txt files with user inputs for defining variables ####
#### reads in a dataset and calculates Z scores relative to given years provided by user ####


def plot_ZvZrepl(inputs_path): 
    
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import PercentFormatter
    from funcs.make_vars import make_vars
    from funcs.make_poly import make_poly
    from funcs.make_poly_mask import make_poly_mask
    from funcs.limit_data_bbox import limit_data_bbox
    from funcs.calc_Z import calc_Z
    from funcs.calc_repl_normal import calc_repl_normal

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

    # Call make_vars function with in_file from .txt file
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
        mask = np.zeros_like(d, dtype=bool)
    
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
    d_mask, ref, Z = calc_Z(d,mask,z0,zf)
    
    # calculate norm_repl arrays for monthly and annual data inputs along with mask for grid cells above user specified threshold
    d90, d_norm_repl, d_ann_norm_repl = calc_repl_normal(d_mask, ref, mask, perc_repl_thresh, 100)

    # create annual data array masked by region polygon     
    d_ann_mask = np.ma.array(d_ann, mask=mask[:len(d_ann)])

    # calculate Z-scores for region data relative to z0 and zf, defined by input file    
    Z = (d_ann_mask - np.mean(d_ann_mask[z0:zf,:,:], axis=0)) / np.std(d_ann_mask[z0:zf,:,:], axis=0)

    # in the same way, calculate Z_repl, simply Z-scores relative to same
    # z0 and zf, but where grid cells >= threshold are replaced by mean 
    #value of cell over reference period
    Z_REPL = (d_ann_norm_repl - np.mean(d_ann_mask[z0:zf,:,:], axis=0)) / np.std(d_ann_mask[z0:zf,:,:], axis=0)
    
    # function to plot distribution of Z vs Z_repl to visualize affects of 
    # grid cells above threshold on annual data
    
    def plot_distribution(Z, Z_REPL, z0, zf):
        
        period = int(zf - z0)

        nmons = len(Z)

        bins = np.linspace(-4, 4, 25)

        plt.figure(figsize=(10,7))
        plt.xlabel('Z')
        plt.ylabel('% Contribution')
        plt.title("{} {} {}".format(dset, dtype, reg_choice))

        #plt.yscale('log')

        ZR1 = Z_REPL[:period,:,:].compressed()
        ZR2 = Z_REPL[nmons-period:,:,:].compressed()

        Z1 = Z[:period,:,:].compressed()
        Z2 = Z[nmons-period:,:,:].compressed()

        plt.hist(ZR1, weights=np.ones(len(ZR1)) / len(ZR1), histtype='step', bins = bins, color = 'w', edgecolor = 'r', linewidth=3, label = '80-93 Z_repl')
        plt.hist(ZR2, weights=np.ones(len(ZR2)) / len(ZR2), histtype='step', bins = bins, color = 'w', edgecolor = 'b', linewidth=3, label = '09-22 Z_repl')

        plt.hist(Z1, weights=np.ones(len(Z1)) / len(Z1), histtype='step', bins = bins, color = 'w', edgecolor = 'r', linestyle = '--', linewidth=1, label = '80-93 Z')
        plt.hist(Z2, weights=np.ones(len(Z2)) / len(Z2), histtype='step', bins = bins, color = 'w', edgecolor = 'b', linestyle = '--', linewidth=1, label = '09-22 Z')

        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.legend(loc='upper left', fontsize=10)
        #plt.xlim(global_min, global_max)
        #plt.gca().invert_yaxis()
        plt.savefig('D:/Research/PRODUCTS/ZvZrepl/{}_{}_ZvZrepl_{}_{}km.png'.format(dset,dtype,reg_choice,grid_res))
        plt.ylim(0, 0.3)
        plt.show()
        
    plot_distribution(Z, Z_REPL, z0, zf)

'''
Z_per = []
# Calculate the period length
period = int(zf - z0)

# Calculate the number of months
nmons = len(Z)

# Define bins for histogram
bins = np.linspace(-3, 3, 25)

# Create a figure
plt.figure(figsize=(13,7))
plt.xlabel('Z')
plt.ylabel('% of dist.')
plt.title('QML Z distribution')

# Iterate over each period
for yr in range(0, 30 -zf):
    start_idx = 0+yr
    end_idx = start_idx+period

    # Extract data for the current period
    Z_period = Z[start_idx:end_idx].compressed()
    Z_per.append(Z_period)

    # Plot histogram for the current period
    hist, bins = np.histogram(Z_period, bins=bins)

    # Plot histogram bars
    plt.hist(Z_period, weights=np.ones(len(Z_period)) / len(Z_period), histtype='step', edgecolor='grey', bins = bins, linewidth=1, alpha=0.3, label='{} - {}'.format(yr+1979,yr+1979+zf))
    
plt.hist(Z[30:].compressed(), weights=np.ones(len(Z_period)) / len(Z_period), histtype='step', edgecolor='b', bins = bins, linewidth=3, label='{} - {}'.format(2009, 2022))

#plt.yscale('log')
plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
plt.legend(loc='upper right', fontsize=7)
# Show plot
plt.show()

##### Calculate Janson Shannon Divergence to compare Z distributions

js_div = []

zf = 15

Z_per = []
# Calculate the period length
period = int(zf - z0)

# Calculate the number of months
nmons = len(Z)

# Define bins for histogram
bins = np.linspace(-3, 3, 25)

hpost, bpost = np.histogram(Z[30:].compressed(), bins=bins)
pde2 = hpost/len(Z[30:].compressed())

# Create a figure
plt.figure(figsize=(13,7))
plt.xlabel('Z')
plt.ylabel('% of dist.')
plt.title('QML Z distribution')

# Iterate over each period
for yr in range(0, 45-zf):
    start_idx = 0+yr
    end_idx = start_idx+period

    # Extract data for the current period
    Z_period = Z[start_idx:end_idx].compressed()
    Z_per.append(Z_period)
    
    plt.figure(figsize=(13,7))

    # Plot histogram for the current period
    hist, bins = np.histogram(Z_period, bins=bins)

    # Plot histogram bars
    plt.hist(Z_period, weights=np.ones(len(Z_period)) / len(Z_period), histtype='step', bins = bins, linewidth=1, label='{} - {}'.format(yr+1979,yr+1979+zf))

    plt.hist(Z[30:].compressed(), weights=np.ones(len(Z_period)) / len(Z_period), histtype='step', edgecolor='b', bins = bins, linewidth=3)

#plt.yscale('log')
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.legend(loc='upper right', fontsize=7)
    
    plt.xlabel('Z')
    plt.ylabel('% of dist.')
    plt.title('QML Z distribution')
    
    plt.gca().set_ylim([0,0.15])
# Show plot
    #plt.savefig('D:/Research/PRODUCTS/z_dists/{}-{}vs{}-{}_DML.png'.format(2009, 2023, yr+1979,yr+1979+zf))
    #plt.close()
    
    #kl_div.append(rel_entr(hpost/len(Z[30:]), hist/len(Z_period)))
    pde1 = hist/len(Z_period)
    #kl_div.append(kl(pde1,pde2))
    js_div.append(distance.jensenshannon(pde1,pde2)**2)
'''
# define directory path where input files are located on device
# looping thru input files for each region , run plot_distribution 
# save .png files to specified location on device
import os    
dir_path = 'D:/Research/SCRIPTS/crich/user_input/Z/ukesm_smb/27km/'
for f_path in os.listdir(dir_path): 
    f = os.path.join(dir_path, f_path)
    plot_ZvZrepl(f)