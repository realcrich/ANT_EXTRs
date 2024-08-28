# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 11:08:13 2024

@author: Collin
"""

# read in netCDF files for CMIP6 participating model 
# (as of 01/24/2024, only CESM2 and UKESM) and interpolate a data product
# (either smb or precip) to a user specified resolution 

def interp_CMIP6_to_polarstereo():
    
    import os
    import numpy as np
    from scipy import interpolate
    import pyproj
    import netCDF4
    
    # set data directory for CMIP6 data
    ddir = 'D:/Research/DATA/CMIP6/'
    
    # prompt user inputs for dataset, datatype, and grid resolution
    dset = input('Which CMIP6 model? (CESM2 or UKESM)?: ')
    dtype = input('smb or pr?: ')
    grid_res = int(input('Resolution of interpolated dataset (km)?: '))
    
    # calculate grid spacing using grid resolution 
    DX = grid_res*1000
    DY = grid_res*1000
    
    def read_in(filepath, variable, grid_res): 

        data = netCDF4.Dataset(filepath)

        # name variables from data file
        time = data.variables['time'][:]
        DATA = data.variables[variable][:,:,:]
        lat = data.variables['lat'][:]
        lon = data.variables['lon'][:]
        
        #-- Grid limits for Bamber DEM
        xmin_bamber=-560.*5e3
        xmax_bamber= 560.*5e3
        ymax_bamber= 560.*5e3
        ymin_bamber=-560.*5e3
        
        #-- create x and y arrays
        x_bamber = np.arange(xmin_bamber,xmax_bamber+DX,DX)
        y_bamber = np.arange(ymax_bamber,ymin_bamber-DY,-DY)

        return DATA, lat, lon, time, x_bamber, y_bamber
    
    def interp_data(data, lon, lat, DX, DY, **kwargs):
        #-- set default keyword arguments for interpolator
        kwargs.setdefault('kx', 1)
        kwargs.setdefault('ky', 1)

        nfiles = data.shape[0]

        #-- coordinate reference systems
        crs1 = pyproj.CRS.from_epsg(4326)
        crs2 = pyproj.CRS.from_epsg(3031)
        transformer = pyproj.Transformer.from_crs(crs1, crs2, always_xy=True)
        direction = pyproj.enums.TransformDirection.INVERSE

        #-- Grid limits for Bamber 1km DEM
        xmin_bamber=-560.*5e3
        xmax_bamber= 560.*5e3
        ymax_bamber= 560.*5e3
        ymin_bamber=-560.*5e3
        #-- create x and y arrays
        x_bamber = np.arange(xmin_bamber,xmax_bamber+DX,DX)
        y_bamber = np.arange(ymax_bamber,ymin_bamber-DY,-DY)
        NX = len(x_bamber)
        NY = len(y_bamber)

        #-- convert x and y axes to grid
        gridx, gridy = np.meshgrid(x_bamber, y_bamber)
        #-- convert from polar stereographic to lat/lon
        gridlon, gridlat = transformer.transform(gridx, gridy, direction=direction)
        #-- allocate for interpolated data
        data_interp = np.ma.zeros((nfiles,NX,NY), fill_value=-9999.0, dtype=np.float64)
        data_interp.mask = np.zeros((nfiles,NX,NY), dtype=bool)

        for mon in range(0,nfiles):
            # create an interpolator for model variable
            s1 = interpolate.RectBivariateSpline(lon, lat, data.data[mon,:,:], **kwargs)
            
            data_interp.data[mon,:,:] = s1.ev(gridlon, gridlat)
            
        return data_interp, x_bamber, y_bamber
    
    def save_ncdf(): 
        
        import netCDF4
        # Assuming data_interp is your interpolated data array
        # loni and lati are the new longitude and latitude arrays

        # Create a new NetCDF file for writing
        output_filepath = 'D:/Research/{}_{}_monthly_{}km_RBS.nc'.format(dset,dtype,grid_res)
        output_nc = netCDF4.Dataset(output_filepath, 'w', format='NETCDF4')

        # Define dimensions in the NetCDF file
        output_nc.createDimension('time', size=None)  # None allows unlimited dimension size
        output_nc.createDimension('lon', size=len(xi))
        output_nc.createDimension('lat', size=len(yi))

        # Create variables in the NetCDF file
        time_var = output_nc.createVariable('time', 'f8', ('time',))
        lon_var = output_nc.createVariable('lon', 'f8', ('lon',))
        lat_var = output_nc.createVariable('lat', 'f8', ('lat',))
        data_var = output_nc.createVariable(dtype, 'f8', ('time', 'lon', 'lat',), fill_value=-9999.0)

        # Write data to NetCDF file
        time_var[:] = np.arange(len(time))  # Assuming time is a 1D array
        lon_var[:] = xi
        lat_var[:] = yi
        data_var[:, :, :] = di  # Assuming di is your interpolated data

        output_nc.close()
    
    # change these input files according to where the CMIP6 data is located on your device
    if dset == 'CESM2': 
        fp1 = os.path.join(ddir,'pr_Amon_CESM2_historical_r1i1p1f1_gn_19790115-20141215_v20190401.nc')
        fp2 = os.path.join(ddir, 'pr_Amon_CESM2_ssp245_r4i1p1f1_gn_20150115-20221215_v20200528.nc')

        fp3 = os.path.join(ddir, 'evspsbl_Amon_CESM2_historical_r1i1p1f1_gn_19790115-20141215_v20190308.nc')
        fp4 = os.path.join(ddir, 'evspsbl_Amon_CESM2_ssp245_r4i1p1f1_gn_20150115-20221215_v20200528.nc')
        
    elif dset == 'UKESM':    
        fp1 = os.path.join(ddir, 'UKESM/pr_Amon_UKESM1-0-LL_historical_r1i1p1f2_gn_19800116-20141216_v20190406.nc')
        fp2 = os.path.join(ddir, 'UKESM/pr_Amon_UKESM1-0-LL_ssp245_r1i1p1f2_gn_20150116-20221216_v20190507.nc')

        fp3 = os.path.join(ddir, 'UKESM/evspsbl_Amon_UKESM1-0-LL_historical_r1i1p1f2_gn_19800116-20141216_v20190627.nc')
        fp4 = os.path.join(ddir, 'UKESM/evspsbl_Amon_UKESM1-0-LL_ssp245_r1i1p1f2_gn_20150116-20221216_v20190715.nc')

    # create seperate precip arrays for the historical data and modelled data
    pr1, lat, lon, time1, x_bamber, y_bamber = read_in(fp1, 'pr', grid_res)
    pr2, lat, lon, time2, x_bamber, y_bamber = read_in(fp2, 'pr', grid_res)
    
    pr = np.concatenate((pr1, pr2), axis=0)
    pr = np.transpose(pr, (0,2,1))
    
    # if precip then interp, if smb calculate by solving dif 
    # between precip and evap, then interp
    if dtype == 'pr':
        
        # concatenate time and data arrays to make full time series across all months
        time = np.concatenate((time1, time2), axis=0)
        
        di, xi, yi = interp_data(pr, lon, lat, DX, DY)
    
    elif dtype == 'smb': 
        
        ev1, lat, lon, time, x_bamber, y_bamber = read_in(fp3, 'evspsbl', grid_res)
        ev2, lat, lon, time, x_bamber, y_bamber = read_in(fp4, 'evspsbl', grid_res)
        
        ev = np.concatenate((ev1, ev2), axis=0)
        ev = np.transpose(ev, (0,2,1))
        
        # concatenate time and data arrays to make full time series across all months
        time = np.concatenate((time1, time2), axis=0)
        
        smb = pr - ev
        
        di, xi, yi = interp_data(smb, lon, lat, DX, DY)
    
    save_ncdf()
    
interp_CMIP6_to_polarstereo()