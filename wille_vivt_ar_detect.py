# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 17:07:16 2024

@author: Collin
"""
import os
import netCDF4 as nc
import xarray as xr
import numpy as np
import datetime
from dateutil.relativedelta import relativedelta
from cytrack.cytrack_functions import get_dates, get_dates_vectors

def run_AR_detection(args_file):
    
    def read_args_from_file(file_path):
        args = {}
        with open(file_path, 'r') as f:
            for line in f:
                key, value = line.strip().split(" = ")
                args[key] = value
        return args

    # Read arguments from the text file
    args = read_args_from_file(args_file)

    # Convert numeric values from strings to their appropriate types
    args['dt_h'] = int(args['dt_h'])
    args['prev_days'] = int(args['prev_days'])

    # Pass arguments to your functions
    idir= args['idir']
    dates, hours = get_dates_vectors(
        year_case_init=args['year_case_init'],
        month_case_init=args['month_case_init'],
        day_case_init=args['day_case_init'],
        hour_case_init=args['hour_case_init'],
        year_case_end=args['year_case_end'],
        month_case_end=args['month_case_end'],
        day_case_end=args['day_case_end'],
        hour_case_end=args['hour_case_end'],
        dt_h=args['dt_h'],
        prev_days=args['prev_days'],
        calendar=args['calendar']
        )

    def get_previous_30_days_files(idir, initial_date):
        date_format = "%Y%m%d"
        hours = ['00', '03', '06', '09', '12', '15', '18', '21']

        start_date = datetime.datetime.strptime(dates[0], date_format)
        end_date = datetime.datetime.strptime(dates[-1], date_format)

        flist = []
        for fp in os.listdir(idir):
            date_hour_str = fp[:-3]  # Assuming the date and hour are in the filename at this position
            try:
                file_datetime = datetime.datetime.strptime(date_hour_str, f"{date_format}_%H")
                if start_date <= file_datetime <= end_date and file_datetime.strftime("%H") in hours:
                    flist.append(os.path.join(idir, fp))
            except ValueError:
                continue  # Skip files that don't match the expected format
        flist = sorted(flist)
        return flist
        
    def calculate_threshold_vivt(flist, lat, lon, lev, g0=9.80665):
        total_vivt = np.zeros((len(flist),len(lat),len(lon)))
        count = 0

        for f in flist:
            ds = xr.open_dataset(f)
            v = ds['V'].values
            qv = ds['QV'].values

            vivt = (1 / g0) * np.trapz(qv * v, x=lev[:, np.newaxis, np.newaxis], axis=0)

            total_vivt[count] = vivt
            count += 1

        thresh = np.nanpercentile(total_vivt,93,axis=0)
        return thresh
    
    def plot_AR_detections(ar_data, vivt, fdate, latg, long, lon, lat): 
        
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import matplotlib.colors as mcolors

        norm = plt.Normalize(0, 1500)
        cmap = mcolors.LinearSegmentedColormap.from_list("", ['white', 'aliceblue', 'cornflowerblue', 'blue', 'limegreen', 'green', 'khaki', 'gold', 'orange', 'orangered', 'firebrick', 'darkorchid'])

        # Create a boolean mask initialized to True
        ar_mask = np.ones_like(vivt, dtype=bool)

        # Set the mask to False at locations corresponding to ar_data
        for j, i in ar_data:
            lat_idx = np.abs(latg[:, 0] - j).argmin()
            lon_idx = np.abs(long[0, :] - i).argmin()
            ar_mask[lat_idx, lon_idx] = False

        ar_mask = ar_mask.astype(int)

        fig = plt.figure(figsize=(12, 12))
        ax = plt.axes(projection=ccrs.SouthPolarStereo())
        ax.set_extent([-180, 180, -90, -40], ccrs.PlateCarree())

        ax.coastlines(resolution="110m", color='k', linewidth=3, zorder=4)
        gl = ax.gridlines(linestyle='--', color='w', draw_labels=True)
        gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 10))
        gl.ylocator = plt.FixedLocator(np.arange(-90, -54, 10))

        lvls = np.linspace(0, 1500, 21)
        cf = plt.contourf(lon, lat, vivt, cmap=cmap, norm=norm, levels=lvls, extend='both', transform=ccrs.PlateCarree(), alpha=0.9)

        c = ax.contour(long, latg, ar_mask, levels=[0.5], colors=['k'], linewidths=2, linestyles='--', transform=ccrs.PlateCarree(), zorder=5)

        cb = plt.colorbar(cf, orientation="vertical", extend='both', fraction=0.046, pad=0.02)
        cb.set_label('vIVT (kg m-1 s-1)', size=8, rotation=90, labelpad=12)
        cb.ax.tick_params(labelsize=10)
        plt.title('AR detections - {}'.format(fdate))
        plt.savefig('D:/Research/DATA/wille_ar_detections/plots/ar_detections_{}.png'.format(fdate))
        plt.close()

    def wille_vivt_ar_detect(idir,dates,hours):

        # create list to store all files from directory
        flist = []
        for fp in os.listdir(idir):
            date_str = fp[:-3]  # Assuming the date is in the filename at this position
            if date_str in [f"{d}_{h}" for d, h in zip(dates, hours)]:
                fname = date_str+'.nc'
                flist.append(os.path.join(idir, fname))
        flist = sorted(flist)

        ex_ds = nc.Dataset(flist[0])

        #t = ex_ds['time'][:]
        lon = ex_ds['lon'][:]
        lat = ex_ds['lat'][:]
        lev = ex_ds.variables['lev'][:] * 100    # Pa --> hPa
        
        lat = lat[::-1]
        lon = lon + 180

        long, latg = np.meshgrid(lon, lat)
        
        coast_mask = np.array(nc.Dataset('D:/Research/DATA/masks/ANT_coast_mask.nc')['ANT_coast_mask'])
        coast_mask = coast_mask[::-1,:].astype(bool)
        
        g0 = 9.80665  # Gravitational acceleration in m/s^2
        #vivt = np.zeros((len(flist)*nt,nlat,nlon))

        # Iterate over daily files
        for f in flist[:1]:
            
            fdate = os.path.basename(f).split('.')[0]
            
            prev_mon_files = get_previous_30_days_files(idir,fdate)
            thresh = calculate_threshold_vivt(prev_mon_files,lat,lon,lev,g0)[::-1,:]
            
            lat_index_list = []
            lon_index_list = []
        #for f in flist[:1]:
            # read the file in and extract daily 3H vIVT
            with nc.Dataset(f) as ds:

                v = ds.variables['V'][:,:,:]
                qv = ds.variables['QV'][:,:,:]

                vivt = (1/g0)*np.trapz(qv * v , x = lev[:,np.newaxis,np.newaxis], axis=0)

                vivt = vivt[::-1, :]
                vivt = np.ma.array(data=vivt,mask=coast_mask)

                res_lon = lon[1] - lon[0]
                res_lat = lat[1] - lat[0]

                origin_lat, origin_lon = latg[0, 0], long[0, 0]
                lat_stepsize = latg[1, 0] - latg[0, 0]
                lon_stepsize = long[0, 1] - long[0, 0]

                ards = np.ma.array(vivt,mask=thresh)
                nards = ards.mask.astype(int)

                indices = np.where(vivt > thresh)
                y_scan = lat[indices[0]]
                x_scan = lon[indices[1]]

                y_splitted_temp = np.split(y_scan, np.where(np.diff(y_scan) < res_lat)[0] + 1)
                x_splitted_temp = np.split(x_scan, np.where(np.diff(y_scan) < res_lat)[0] + 1)

                y_longest = max(y_splitted_temp, key=len)
                x_longest = max(x_splitted_temp, key=len)

                x_reverse = []
                y_reverse = []

                try:
                    if y_longest.max() - y_longest.min() > 20:
                        reverse_grid = np.arange(min(x_longest), max(x_longest) + 0.5, res_lon)
                        for i in reverse_grid:
                            x_index_reverse = np.where(x_longest == i)
                            x_reverse = np.concatenate((x_reverse, x_longest[x_index_reverse]))
                            y_reverse = np.concatenate((y_reverse, y_longest[x_index_reverse]))
                except ValueError:
                    pass

                x_splitted = np.split(x_reverse, np.where(np.diff(x_reverse) > 20)[0] + 1)
                y_splitted = np.split(y_reverse, np.where(np.diff(x_reverse) > 20)[0] + 1)

                try:
                    if x_splitted[0][0] + 360 - x_splitted[-1][-1] < 20:
                        x_splitted[-1] = np.concatenate((x_splitted[-1], x_splitted[0]))
                        x_splitted = np.delete(x_splitted, 0, 0)
                        y_splitted[-1] = np.concatenate((y_splitted[-1], y_splitted[0]))
                        y_splitted = np.delete(y_splitted, 0, 0)
                except IndexError:
                    pass

                x_shape = []
                y_shape = []

                for i in range(len(y_splitted)):
                    x_reverse2 = []
                    y_reverse2 = []
                    x_final = []
                    y_final = []

                    try:
                        reverse_grid2 = np.arange(max(y_splitted[i]), min(y_splitted[i]) - 0.5, res_lat)
                        for j in reverse_grid2:
                            y_index_reverse2 = np.where(y_splitted[i] == j)
                            x_reverse2 = np.concatenate((x_reverse2, x_splitted[i][y_index_reverse2]))
                            y_reverse2 = np.concatenate((y_reverse2, y_splitted[i][y_index_reverse2]))
                    except ValueError:
                        pass

                    try:
                        y_splitted_final = np.split(y_reverse2, np.where(np.diff(y_reverse2) < res_lat)[0] + 1)
                        x_splitted_final = np.split(x_reverse2, np.where(np.diff(y_reverse2) < res_lat)[0] + 1)

                        for z in range(len(y_splitted_final)):
                            if y_splitted_final[z].max() - y_splitted_final[z].min() > 20:
                                y_final = np.concatenate((y_final, y_splitted_final[z]))
                                x_final = np.concatenate((x_final, x_splitted_final[z]))

                        x_shape = np.concatenate((x_shape, x_final))
                        y_shape = np.concatenate((y_shape, y_final))
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

                    # Save AR detection data to a file in the specified format
                ar_data = []
                #ar_lat=lat[::-1]
                ar_lon=lon-180
                for lat_idx, lon_idx in zip(lat_index_list, lon_index_list):
                    if len(lat_idx) > 0:
                        for j,i in zip(lat_idx, lon_idx):
                            # Check if the point (j, i) is outside the coast_mask
                            if coast_mask[j, i] == False:
                                ar_data.append([lat[j], lon[i]])
                            else: 
                                pass
                                
                ar_data = np.array(ar_data)
                np.savetxt(f"D:/Research/DATA/wille_ar_detections/trajectories/ar_trajectories_{fdate}.dat", ar_data, fmt='%d', header="Latitude_Index Longitude_Index")
                plot_AR_detections(ar_data, vivt, fdate, latg, long, lon, lat)
        
    wille_vivt_ar_detect(idir,dates,hours)
    #return lat_index_list, lon_index_list