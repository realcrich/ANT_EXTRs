# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 13:19:34 2023

@author: Collin
"""

import cartopy.crs as ccrs
import numpy as np
import h5py
import glob
import os

def make_ATL11_vars(ATL11_dir, dh):
    
    atl11_dir = 'D:/Research/DATA/ICESAT2/ATL11/'
    
    def read_field(pt, var): 
        data = np.array(pt[var])
        bad = (data==pt[var].attrs['_FillValue'])
        data[bad]=np.NaN
        return data
        
    def extract_vars(file): 
        with h5py.File(file, 'r') as h5f:
            for pair in ['pt1','pt2','pt3']:
                lat = np.array(h5f[pair]['latitude'])
                lon = np.array(h5f[pair]['longitude'])
                
                dt = np.array(h5f[pair]['delta_time'])
                
                h_corr = np.array(h5f[pair]['h_corr'])
                h_corr_sig = np.array(h5f[pair]['h_corr_sigma'])
                h_corr_sig_s = read_field(h5f[pair], 'h_corr_sigma_systematic')
                qual = np.array(h5f[pair]['quality_summary'])            
            for col in range(h_corr.shape[1]):
                h_corr[qual==6]=np.NaN
                
        h_corr[h_corr == h_corr.max()] = np.NaN
        h_corr_sig[h_corr_sig == h_corr_sig.max()] = np.NaN
        h_corr_sig_s[h_corr_sig_s == h_corr_sig_s.max()] = np.NaN
        
        return lat, lon, h_corr, np.sqrt(h_corr_sig**2+h_corr_sig_s**2), dt

    files=sorted(glob.glob(os.path.join(ATL11_dir, 'ATL11*.h5')))
   
    lon=[]
    lat=[]
    h_corr=[]
    sigma_h=[]
    dt = []
    ncycs = 0
    for file in files:
        lats, lons, hh, ss, ts = extract_vars(file)
        lon += [lons]
        lat += [lats]
        h_corr += [hh]
        sigma_h += [ss]
        dt += [ts]
        ncycs+=1
        
        print(file)
    
    lon=np.concatenate(lon)
    lat=np.concatenate(lat)
    h_corr=np.concatenate(h_corr, axis=0)
    sigma_h=np.concatenate(sigma_h, axis=0)
    dt=np.concatenate(dt, axis=0)

    projection = ccrs.SouthPolarStereo()

    return ncycs, lon, lat, h_corr, sigma_h, dt, projection

    if dh == 'y':
        from funcs.read_ATL14 import read_ATL14
        
        lon14, lat14, h_ref, h_ref_sigma, projection = read_ATL14()
        
        dh_atl11 = np.zeros_like(h_corr)

        for per in range(ncycs): 
            dh_atl11[:,per] = h_corr[:,per] - h_ref
        return dh_atl11, lon, lat, dt
    else:
        return h_corr, lon, lat, dt

    lon, lat, h, h_sig, dt, proj = make_ATL11_vars('D:\\Research\\DATA\\ICESAT2\\ATL11\\','y')
