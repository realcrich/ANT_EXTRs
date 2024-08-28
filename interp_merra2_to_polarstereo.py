# -*- coding: utf-8 -*-
"""
Created on Wed Jan 17 12:49:35 2024

@author: Collin
"""

# read in netCDF file for MERRA2 reanalysis and interpolate a data product
# (either smb or precip) to a user specified resolution 

def interp_merra2_to_polarstereo(): 
    
    import os
    import netCDF4
    import numpy as np
    from scipy import interpolate
    import pyproj
    from funcs.make_vars import make_vars
    
    # set data directory for CMIP6 data
    ddir = 'D:/Research/DATA/MERRA2/'

    def interp_data(data, glon, glat, DX, DY, **kwargs):
        
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

        for mon in range(0,nfiles):
            # create an interpolator for model variable
            s1 = interpolate.RectBivariateSpline(glon, glat, data.data[mon,:,:], **kwargs)
            data_interp.data[mon,:,:] = s1.ev(gridlon, gridlat)

        return data_interp, x_bamber, y_bamber
    
    def save_ncdf(): 
        # Assuming data_interp is your interpolated data array
        # loni and lati are the new longitude and latitude arrays

        # Create a new NetCDF file for writing
        output_filepath = 'D:/Research/MERRA2_{}_monthly_{}km_RBS.nc'.format(dtype,grid_res)
        output_nc = netCDF4.Dataset(output_filepath, 'w', format='NETCDF4')

        # Define dimensions in the NetCDF file
        output_nc.createDimension('time', size=None)  # None allows unlimited dimension size
        output_nc.createDimension('lon', size=len(loni))
        output_nc.createDimension('lat', size=len(lati))

        # Create variables in the NetCDF file
        time_var = output_nc.createVariable('time', 'f8', ('time',))
        lon_var = output_nc.createVariable('lon', 'f8', ('lon',))
        lat_var = output_nc.createVariable('lat', 'f8', ('lat',))
        data_var = output_nc.createVariable(dtype, 'f8', ('time', 'lon', 'lat',), fill_value=-9999.0)

        # Write data to NetCDF file
        time_var[:] = np.arange(len(time))  # Assuming time is a 1D array
        lon_var[:] = loni
        lat_var[:] = lati
        data_var[:, :, :] = di  # Assuming di is your interpolated data

        output_nc.close()
    
    # prompt user for inputs for MERRA2 data filepath, datatype, and grid res
    fpath = input('Filepath to dataset?: ')
    f = os.path.join(ddir,fpath)
    dtype = input('smb or precip?: ')
    grid_res = int(input('Resolution of interpolated dataset (km)?: '))
    
    dx = grid_res*1000
    dy = grid_res*1000
    
    d, d_ann, time, lon, lat, x_grid, y_grid = make_vars(f, dtype, grid_res)
    
    di, loni, lati = interp_data(d, lon, lat, dx, dy)
    
    save_ncdf()
    
    return di, loni, lati

interp_merra2_to_polarstereo()