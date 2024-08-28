# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 10:54:44 2024

@author: Collin
"""

# in_file = specific file path to dataset
# data_type = what variable is being read in ('SMB', 'precip', ..., etc)
# grid_res = (x, y) resolution (in km) of Bamber Polar Stereographic grid (used to make x,y meshgrid of coordinates)

import netCDF4 as nc
import numpy as np

def make_vars(in_file, data_type, grid_res): 

    # read in dataset using in-file
    data = nc.Dataset(in_file)

    # name variables from data file
    time = data.variables['time'][:]
    DATA = data.variables[data_type][:,:,:]
    lat = data.variables['lat'][:]
    lon = data.variables['lon'][:]
    
    num_months = DATA.shape[0]
    num_years = num_months // 12
    
    # make a truncated data array that is an integer multiple of 12, so only complete years are included in DATA_ANN
    truncated_data = DATA[:num_years * 12]

    DATA_ANN = np.ma.sum(truncated_data.reshape(-1, 12, DATA.shape[1], DATA.shape[2]), axis=1)

    xmin_bamber=-560.*5e3
    xmax_bamber= 560.*5e3
    ymax_bamber= 560.*5e3
    ymin_bamber=-560.*5e3
    
    dx = grid_res * 1000
    dy = grid_res * 1000

    x_bamber = np.arange(xmin_bamber,xmax_bamber+dx,dx)
    y_bamber = np.arange(ymax_bamber,ymin_bamber-dy,-dy)

    x_grid, y_grid = np.meshgrid(x_bamber, y_bamber)

    return DATA, DATA_ANN, time, lon, lat, x_grid, y_grid