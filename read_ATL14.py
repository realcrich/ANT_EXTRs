# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 13:07:15 2023

@author: Collin
"""

import cartopy.crs as ccrs
import numpy as np
import h5py
import glob
import os

def read_ATL14(): 
    def read_field(pt, var): 
        data = np.array(pt[var])
        bad = (data==pt[var].attrs['_FillValue'])
        data[bad]=np.NaN
        return data

    def read_ATL14(file): 
        with h5py.File(file, 'r') as h5f:
            for pair in ['pt1','pt2','pt3']:
                lat = np.array(h5f[pair]['latitude'])
                lon = np.array(h5f[pair]['longitude'])
                dem_h = np.array(h5f[pair]['ref_surf']['dem_h'])
                dem_h_sigma = np.array(h5f[pair]['ref_surf']['dem_h_sigma'])
                
        dem_h[dem_h == dem_h.max()] = np.NaN
        dem_h_sigma[dem_h_sigma == dem_h_sigma.max()] = np.NaN
                      
        return lat, lon, dem_h, dem_h_sigma

    def make_vars(ATL14_dir):
    	#ATL11_dir="C:\\Users\\Collin\\Documents\\ESS21-22\\Research\\"
    	files=sorted(glob.glob(os.path.join(ATL14_dir, 'ATL14*.h5')))

    	lon=[]
    	lat=[]
    	dem_h=[]
    	dem_h_sigma=[]
    	for file in files:
        		lats, lons, hh, hhs = read_ATL14(file)
        		lon += [lons]
        		lat += [lats]
        		dem_h += [hh]
        		dem_h_sigma += [hhs]
        
    	lon=np.concatenate(lon)
    	lat=np.concatenate(lat)
    	h_ref=np.concatenate(dem_h, axis=0)
    	h_ref_sigma=np.concatenate(dem_h_sigma, axis=0)

    	projection = ccrs.SouthPolarStereo()

    	return lon, lat, h_ref, h_ref_sigma, projection

    lon, lat, h_ref, h_ref_sigma, projection = make_vars('D:/Research/DATA/ICESAT2/ATL14/')
        
    return lon, lat, h_ref, h_ref_sigma, projection

#lons, lats, h_ref, h_ref_sigma, proj = make_vars('D:/Research/DATA/ICESAT2/ATL14/')
#lons, lats, h_ref, h_ref_sigma = make_vars('D:/Research/DATA/ICESAT2/ATL14/')

'''
dh_atl11 = np.zeros_like(h)

for per in range(17): 
    dh_atl11[:,per] = h[:,per] - h_ref
    

dt[dt == np.max(dt)] = np.NaN

from astropy.time import Time

t18 = Time('2018-01-01T00:00:00')

t20 = Time('2020-01-01T00:00:00')

t80=Time('1980-01-06T00:00:00')

t18 = t18.to_value('gps')

t20 = t20.to_value('gps')

t80 = t80.to_value('gps')

t_dif = t20 - t18

t_anom = dt - t_dif
'''