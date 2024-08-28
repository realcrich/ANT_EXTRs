# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 16:39:08 2024

@author: Collin
"""

def interp_racmo2_to_polarstereo(): 
    
    import os
    import netCDF4
    import numpy as np
    from scipy import interpolate
    import pyproj
    
    # set data directory for CMIP6 data
    #ddir = 'D:/Research/DATA/RACMO/RACMO_Updated_may_2022.dir/RACMO_Updated_may_2022.dir/'
    ddir = 'D:/Research/DATA/RACMO/197901-202312/'
    
    def read_in_racmo(filepath, variable): 

        data = netCDF4.Dataset(filepath)
                    
        # name variables from data file
        time = data.variables['time'][:]
        DATA = data.variables[variable][:,:,:,:].squeeze()
        lat = data.variables['lat'][:,:]
        lon = data.variables['lon'][:,:]
        rlat = data.variables['rlat'][:]
        rlon = data.variables['rlon'][:]

                
        return DATA, lat, lon, time, rlat, rlon

    def interp_data(data, glon, glat, DX, DY, **kwargs):
        
        #-- set default keyword arguments for interpolator
        kwargs.setdefault('kx', 1)
        kwargs.setdefault('ky', 1)

        nfiles = data.shape[0]

        #-- coordinate reference systems
        crs1 = pyproj.CRS.from_string('-m 57.295779506 +proj=ob_tran +o_proj=latlon +o_lat_p=-180.0 +lon_0=10.0') 
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
            s1 = interpolate.RectBivariateSpline(glon, glat, data[mon,:,:].T, **kwargs)

            data_interp.data[mon,:,:] = s1.ev(gridlon, gridlat)

        return data_interp, x_bamber, y_bamber
    
    def save_ncdf(): 
        # Assuming data_interp is your interpolated data array
        # loni and lati are the new longitude and latitude arrays

        # Create a new NetCDF file for writing
        output_filepath = 'D:/Research/RACMO2.4_196001-202212_{}_monthly_{}km_RBS.nc'.format(dtype,grid_res)
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
    
    # prompt user inputs for dataset, datatype, and grid resolution
    fpath = input('Filepath to dataset?: ')
    f = os.path.join(ddir,fpath)
    dtype = input('smb or precip?: ')
    grid_res = int(input('Resolution of interpolated dataset (km)?: '))
    
    dx = grid_res*1000
    dy = grid_res*1000
    
    d, lat, lon, time, rlat, rlon = read_in_racmo(f, dtype)
    
    di, loni, lati = interp_data(d, rlon, rlat, dx, dy)
    
    save_ncdf()
    
    return di, loni, lati