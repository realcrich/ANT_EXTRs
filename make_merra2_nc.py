# -*- coding: utf-8 -*-
"""
Created on Fri Jan 12 11:08:31 2024

@author: Collin
"""

# reads in all monthly merra2 .nc files acquired via Tyler Sutterley's gesdisc_merra_sync.py script
# and calculates a user specified product (smb or precip), and makes a single .nc file with all
# months included to use in interpolation routine, interp_merra2_to_polarstereo.py  

def make_merra2_nc():
    
    import os
    import numpy as np
    import netCDF4 as nc
    
    # set merra2 data directory
    merra2_dir = 'D:/Research/DATA/MERRA2' 
    
    # prompt user inputs for desired product and HEM, hemisphere, for either ANT or GRN
    dtype = input('Product to be derived (precip or smb)?: ')
    HEM = input('North or South Pole (N/S)?: ')
    
    # set empty lists for file names  
    int_files = []
    glc_files = []
    products = []
    nmons = 0
    
    dint_dir = os.path.join(merra2_dir,'M2TMNXINT.5.12.4')
    for root, dirs, files_in_dir in os.walk(dint_dir):
        for file in files_in_dir:
            if file.endswith('.nc4'):
                file_path = os.path.join(root, file)
                int_files.append(file_path)
                nmons += 1
                
    if dtype == 'precip':
        products.append(('PRECCU', 'PRECCLS', 'PRECSN'))
        
    elif dtype == 'smb':
        dglc_dir = os.path.join(merra2_dir, 'M2TMNXGLC.5.12.4')
        for root, dirs, files_in_dir in os.walk(dglc_dir):
            for file in files_in_dir:
                if file.endswith('.nc4'):
                    file_path = os.path.join(root, file)
                    glc_files.append(file_path)
        products.append(('PRECCU', 'PRECCLS', 'PRECSN', 'EVAP', 'RUNOFF', 'WESNSC'))
    
    # slice the first month of data to make vars for the length of lon and lat arrays     
    ex_ds = nc.Dataset(int_files[0])
        
    if HEM == 'S':
        lon = ex_ds.variables['lon'][:]
        lat = ex_ds.variables['lat'][:61]
    
        nlon = len(lon)
        nlat = len(lat)
    elif HEM == 'N':
        lon = ex_ds.variables['lon'][:]
        lat = ex_ds.variables['lat'][:]
    
        nlon = len(lon)
        nlat = len(lat)
        
    # create arrays for all monthly data and time array for each month
    d = np.ma.zeros((nmons, nlon, nlat), fill_value=-9999.0, dtype=np.float64)
    t = np.zeros((nmons))
    dpm = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]  # days per minute
    spd = 24.*60.*60.    # seconds per day
    
    for mon in range(nmons): 
        
        ds_int = nc.Dataset(int_files[mon])
        moy = int(int_files[mon][-6:-4]) - 1
        dfloat = dpm[moy]*spd / 365.25
        
        PRECCU = ds_int.variables['PRECCU'][:,:nlat,:].squeeze().T * dpm[moy]*spd
        PRECLS = ds_int.variables['PRECLS'][:,:nlat,:].squeeze().T * dpm[moy]*spd
        PRECSN = ds_int.variables['PRECSN'][:,:nlat,:].squeeze().T * dpm[moy]*spd

        if dtype == 'precip':
            d[mon,:,:] = PRECCU + PRECLS + PRECSN 
            t[mon] = dfloat
            
        elif dtype == 'smb': 
            ds_glc = nc.Dataset(glc_files[mon])
            EVAP = ds_int.variables['EVAP'][:,:nlat,:].squeeze().T * dpm[moy]*spd
            RUNOFF = ds_glc.variables['RUNOFF'][:,:nlat,:].squeeze().T * dpm[moy]*spd
            WESNSC = ds_glc.variables['WESNSC'][:,:nlat,:].squeeze().T * dpm[moy]*spd
            
            d[mon,:,:] = PRECCU + PRECLS + PRECSN - EVAP - RUNOFF - WESNSC
            t[mon] = dfloat
            
    def save_ncdf(): 
        # Assuming data_interp is your interpolated data array
        # loni and lati are the new longitude and latitude arrays

        # Create a new NetCDF file for writing
        output_filepath = 'D:/Research/MERRA2_{}_monthly.nc'.format(dtype)
        output_nc = nc.Dataset(output_filepath, 'w', format='NETCDF4')

        # Define dimensions in the NetCDF file
        output_nc.createDimension('time', size=len(t))  # None allows unlimited dimension size
        output_nc.createDimension('lon', size=len(lon))
        output_nc.createDimension('lat', size=len(lat))

        # Create variables in the NetCDF file
        time_var = output_nc.createVariable('time', 'f8', ('time',))
        lon_var = output_nc.createVariable('lon', 'f8', ('lon',))
        lat_var = output_nc.createVariable('lat', 'f8', ('lat',))
        data_var = output_nc.createVariable(dtype, 'f8', ('time', 'lon', 'lat',), fill_value=-9999.0)

        # Write data to NetCDF file
        time_var[:] = t  # Assuming time is a 1D array
        lon_var[:] = lon
        lat_var[:] = lat
        data_var[:, :, :] = d  # Assuming di is your interpolated data

        output_nc.close()
        
    save_ncdf()
        
make_merra2_nc()