# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 12:04:37 2024

@author: Collin
"""

import numpy as np
from funcs.map_ll import map_ll
from shapely.ops import cascaded_union
from funcs.make_poly import make_poly
import os
import netCDF4 as nc
from cytrack.cytrack_functions import *


def wille_ivt_ar_detect(ddir,month,year,outdir,regs):
    
    # create list to store all files from directory
    flist = []
    for fp in os.listdir(ddir):
        if np.logical_and(str(fp[-7:-5]) == month, str(fp[-11:-7]) == year):
            flist.append(os.path.join(ddir,fp))
    flist = sorted(flist)

    # Define months and corresponding number of days
    mons = ['01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12']
    days = ['31', '28', '31', '30', '31', '30', '31', '31', '30', '31', '30', '31']
    years = ['2002', '2003', '2004', '2005', '2006', '2007', '2008', '2009', '2010', \
         '2011', '2012', '2013','2014', '2015', '2016', '2017', '2018', '2019', \
             '2020', '2021', '2022', '2023']

    #regs = ['aap', 'apb', 'ka']
    reg_choices = ['aap', 'apb', 'bc', 'ccp', 'cpd', 'ddp', 'dpe', 'eep', 'epf', 'fg', 'gh', 'hhp', 'hpi', 'iipp', 'ippj', 'jjpp', 'jppk', 'ka']

    # make empty list to store polygon for each basin
    polys = []

    # read in .shp file and exract region poly, append to list
    for reg in regs:
        _,_,_,polyreg,_,_,_,_ = make_poly(reg,'n')
        polys.append(polyreg)

    # combine all reg polys into single poly for all Antarctica
    polyc = cascaded_union(polys)

    ex_ds = nc.Dataset(flist[0])

    t = ex_ds['time'][:]
    lon = ex_ds['lon'][:]
    lat = ex_ds['lat'][:]
    lev = ex_ds.variables['lev'][:] * 100    # Pa --> hPa

    nt = len(ex_ds['time'][:])
    nlev = len(ex_ds['lev'][:]) 
    nlon = len(ex_ds['lon'][:])
    nlat = len(ex_ds['lat'][:])
    
    # create meshgrid of lon, lat coords to use with polygon
    long, latg = np.meshgrid(lon, lat)
    
    # convert to polar stereographic projection
    x,y = map_ll(long, latg, HEM='S', FORMAT='tuple')

    # initialize mask for ANT, mask any grid cell that is within grounded ice sheet
    points = geometry.MultiPoint(np.column_stack((x.flatten(), y.flatten())))
    coast_mask = np.array([polyc.contains(point) for point in points]).reshape(len(lat), len(lon))[::-1,:]

    g0 = 9.80665  # Gravitational acceleration in m/s^2
    vivt = np.zeros((len(flist)*nt,nlat,nlon))

    # Iterate over daily files
    for f in flist:

        idx = 0

            # read the file in and extract daily 3H vIVT
        with nc.Dataset(f) as ds:
            
            for tidx in range(nt): 

                v = ds.variables['V'][tidx,:,:,:]
                qv = ds.variables['QV'][tidx,:,:,:]

                vIVT = (1/g0)*np.trapz(qv * v , x = lev[:,np.newaxis,np.newaxis], axis=0)

                vivt[idx] = vIVT

                idx += 1
        
    vivt = vivt[:,::-1,:]
    lat = lat[::-1]
    lon = lon+180

    # create meshgrid of lon, lat coords to use with polygon
    long, latg = np.meshgrid(lon, lat)

    res_lon = lon[1] - lon[0]
    res_lat = lat[1] - lat[0]

    origin_lat,origin_lon = latg[0,0],long[0,0]
    lat_stepsize = latg[1,0] - latg[0,0]
    lon_stepsize = long[0,1] - long[0,0]
    
    lat_index_list = []
    lon_index_list = []
    landfall_idx = []
    thresh = np.zeros_like(vivt)    # kg m -1 s-1
    thresh[:,:,:] = 100

    indices = np.where(vivt > thresh)
    time = indices[0]
    y_scan = lat[indices[1]]
    x_scan = lon[indices[2]]

    landfall_idx = []
    timesteps = np.arange(0,len(vivt))

    ar_mask = np.zeros_like(coast_mask)
    ar_mask = np.tile(ar_mask, len(timesteps)).reshape(len(timesteps),len(lat),len(lon))
    landfall_mask = np.zeros_like(ar_mask)

    for timestep in timesteps :
        #print(timestep)
        timestep_idx = np.where(timestep == time)[0]
        y_idx = y_scan[timestep_idx]
        x_idx = x_scan[timestep_idx]


        y_splitted_temp = np.split(y_idx, np.where(np.diff(y_idx) < res_lat)[0] +1)
        x_splitted_temp = np.split(x_idx, np.where(np.diff(y_idx) < res_lat)[0] +1)


        y_longest = max(y_splitted_temp, key=len)
        x_longest = max(x_splitted_temp, key=len)


        x_reverse = []
        y_reverse = []


        try:

            if y_longest.max() - y_longest.min() > 20 :
                reverse_grid = np.arange(min(x_longest),max(x_longest)+0.5,res_lon)
                for i in reverse_grid:
                    x_index_reverse = np.where(x_longest == i)
                    x_reverse = np.concatenate((x_reverse,x_longest[x_index_reverse]))
                    y_reverse = np.concatenate((y_reverse,y_longest[x_index_reverse]))
        except ValueError:
            pass

        x_splitted = np.split(x_reverse, np.where(np.diff(x_reverse) > 20)[0] +1) #This is where the error occurs. If x_reverse contains 358 and 2 for instance, than the shape b
        y_splitted = np.split(y_reverse, np.where(np.diff(x_reverse) > 20)[0] +1) # If x_reverse contains 358 and 2 for instance, than the shape is incorrectly split




        try:
            if x_splitted[0][0]+360 - x_splitted[-1][-1] < 20:
                x_splitted[-1] = np.concatenate((x_splitted[-1],x_splitted[0]))
                x_splitted = np.delete(x_splitted,0,0)
                y_splitted[-1] = np.concatenate((y_splitted[-1],y_splitted[0]))
                y_splitted = np.delete(y_splitted,0,0)


        except IndexError:
            pass
        x_shape = []
        y_shape = []

        x_shape_landfall = []
        y_shape_landfall = []

        for i in range(0,len(y_splitted)):
            x_reverse2 = []
            y_reverse2 = []

            x_final = []
            y_final = []


            try:	
                reverse_grid2 = np.arange(max(y_splitted[i]),min(y_splitted[i])-0.5,res_lat)
                for j in reverse_grid2:
                    y_index_reverse2 = np.where(y_splitted[i] == j)
                    x_reverse2 = np.concatenate((x_reverse2,x_splitted[i][y_index_reverse2]))
                    y_reverse2 = np.concatenate((y_reverse2,y_splitted[i][y_index_reverse2]))

            except ValueError:
                pass
            try: 
                y_splitted_final = np.split(y_reverse2, np.where(np.diff(y_reverse2) < res_lat)[0] +1)
                x_splitted_final = np.split(x_reverse2, np.where(np.diff(y_reverse2) < res_lat)[0] +1)

                for z in range(0, len(y_splitted_final)):
                    if y_splitted_final[z].max() - y_splitted_final[z].min() > 20:
                        y_final = np.concatenate((y_final, y_splitted_final[z]))
                        x_final = np.concatenate((x_final, x_splitted_final[z]))

                x_shape = np.concatenate((x_shape, x_final))
                y_shape = np.concatenate((y_shape,y_final))

            except ValueError:
                pass
        lat_index = (y_shape - origin_lat) / lat_stepsize
        lon_index = (x_shape - origin_lon) / lon_stepsize
        lat_index = lat_index.astype(int)
        lon_index = lon_index.astype(int)

        if len(lat_index) > 0:
            lat_index_list.append(lat_index.data)
        else:
            lat_index_list.append(lat_index)

        if len(lon_index) > 0:
            lon_index_list.append(lon_index.data)
        else:
            lon_index_list.append(lon_index)

        if len(lat_index_list[timestep]) > 0:
            ar_mask[timestep,lat_index_list[timestep],lon_index_list[timestep]] = True
            landfall_mask[timestep] = np.logical_and(ar_mask[timestep],coast_mask)
            if np.any(landfall_mask[timestep] == True):
                landfall_idx.append(((timestep+1)/8, np.where(landfall_mask[timestep] == True)))
            else:
                landfall_idx.append([])
        else:
            landfall_idx.append([])
            
    return lat_index_list, lon_index_list, landfall_idx, ar_mask, landfall_mask