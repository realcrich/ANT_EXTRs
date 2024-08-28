# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:32:35 2024

@author: Collin
"""

import os
import netCDF4 as nc
import xarray as xr
import numpy as np
import datetime
from datetime import timedelta
import pandas as pd
from dateutil.relativedelta import relativedelta
from cytrack.cytrack_functions import get_dates, get_dates_vectors

def read_args_from_file(file_path):
    
    """
    Read arguments from parameters input file : 'wille_ar_detections_inputs.txt'

    Parameters:
    file_path (dtype: str): name of file containing input parameters for run_ARDetect.py

    Returns:
    args (dtype: list): arguments input parameters.
    """
    
    args = {}
    with open(file_path, 'r') as f:
        for line in f:
            key, value = line.strip().split(" = ")
            args[key] = value
    return args

# Define function to find AR landfall points
def find_ar_landfall_points(ar_mask,lat,lon):
    
    """
    find cells above the AR threshold that make landfall according to ar_mask

    Parameters:
    ar_mask (dtype: masked array): 2D numpy array for defining which points meet AR criteria
    lat (dtype: numpy.ndarray) : 1D array of latitude values
    lon (dtype: numpy.ndarray) : 1D array of longitude values

    Returns:
    landfall_points (dtype: list): longitudinal points where AR meets the coastline.
    """
    
    landfall_points = []
    for j in range(1, ar_mask.shape[0] - 1):
        for i in range(1, ar_mask.shape[1] - 1):
            if ar_mask.data[j, i] == 0:
                # Check if any of the 8 neighbors are True in the mask
                neighbors = [
                    ar_mask.mask[j-1, i], ar_mask.mask[j+1, i], ar_mask.mask[j, i-1], ar_mask.mask[j, i+1],
                    ar_mask.mask[j-1, i-1], ar_mask.mask[j-1, i+1], ar_mask.mask[j+1, i-1], ar_mask.mask[j+1, i+1]
                ]
                if any(neighbors):
                    landfall_points.append((lat[j], lon[i]))
    return landfall_points

# Function to find contiguous groups of longitudes
def find_contiguous_groups(lons):
    
    """
    find continous groups of points where the AR makes landfall along the coast ; accepts 
    input from find_ar_landfall_points(). 

    Parameters:
    lons (dtype: list): longitudes where AR makes landfall

    Returns:
    groups (dtype: list): list of lists, where each is a continuous segment of longitudes.
    """
    
    if len(lons)>0:
        sorted_lons = sorted(lons)
        groups = []
        group = [sorted_lons[0]]

        for lon in sorted_lons[1:]:
            if lon - group[-1] <= 1:
                group.append(lon)
            else:
                groups.append(group)
                group = [lon]
        groups.append(group)
    else:
        groups =[]
    return groups

def get_previous_30_days_files(idir, dates, initial_date):
    
    """
    find datetime objects for 30 days prior to initial date in dates. 

    Parameters:
    idir (dtype: list): directory where data files are located
    dates (dtype: list): list of datetime objects as specified by arguments input file
    initial_date (dtype: datetime object): first YYYYMMDD-HH in dates

    Returns:
    flist (dtype: list): list of all files, defined as valid via dates.
    """
    
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
    
    """
    Calculate a threshold vIVT mask, to represent values that need be exceeded to be considered an AR.
    Threshold values for each point are calculated using the mean over the previous 30 days.

    Parameters:
    flist (dtype: list): list of files to be read in to calculate the threshold vIVT 
    lat (dtype: numpy.ndarray) : 1D array of latitude values
    lon (dtype: numpy.ndarray) : 1D array of longitude values
    lev (dtype: numpy.ndarray) : 1d array of pressure levels in Pa
    g0 (dtype: float) : gravitational acceleration constant

    Returns:
    thresh (dtype: array): 2D array of vIVT values above the 93rd percentile relative to the previous 30 days
    """
    
    total_vivt = np.zeros((len(flist),len(lat),len(lon)))
    count = 0

    for f in flist:
        ds = xr.open_dataset(f)
        v = ds['V'].values
        #u = ds['U'].values
        qv = ds['QV'].values

        qvv = np.trapz(qv * v, x=lev[:, np.newaxis, np.newaxis], axis=0)
        #qvu = np.trapz(qv * u, x=lev[:, np.newaxis, np.newaxis], axis=0)
        #vivt = (1 / g0) * np.sqrt(qvu**2+qvv*2)
        vivt = (1/g0) * qvv

        total_vivt[count] = vivt
        count += 1

    thresh = np.nanpercentile(total_vivt,93,axis=0)[::-1,:]
    return thresh

def calc_dx_dy(lat, lon):
    
    """
    Convert latitude and longitude grid spacings to meters.

    Parameters:
    lat (numpy.ndarray): 1D array of latitudes.
    lon (numpy.ndarray): 1D array of longitudes.

    Returns:
    dlat_m, dlon_m (dtype: tuple): Grid spacings in meters (dlat_m, dlon_m).
    """
    
    # Radius of the Earth in meters
    R = 6371000

    # Convert latitude to meters
    dlat = np.gradient(lat)
    dlat_m = dlat * (np.pi / 180) * R

    # Convert longitude to meters
    dlon = np.gradient(lon)
    lat_rad = np.radians(lat)
    dlon_m = np.zeros((len(lat), len(dlon)))

    for i in range(len(lat)):
        dlon_m[i, :] = dlon * (np.pi / 180) * R * np.cos(lat_rad[i])

    return dlat_m, dlon_m

def calculate_diabatic_heating(f,lat,lon,lev):
    
    """
    Calculate diabatic heating rate.

    Parameters:
    f (dtype: str): filename for input data file.
    lat (dtype: numpy.ndarray): 1D array of latitudes.
    lon (dtype: numpy.ndarray): 1D array of longitudes.
    lev (dtype: numpy.ndarray) : 1d array of pressure levels in Pa

    Returns:
    theta, Q, viQ, Qkat (dtype: tuple): theta (potential temperature), Q (diabatic heating rate), 
                                        viQ (vertically integrated diabatic heating rate), 
                                        Qkat (diabatic heating from only x-component).
    """
    
    R = 287.05  # Specific gas constant for dry air in J/(kg*K)
    Cp = 1005 # Specific heat at constant pressure for air J/(kg*K)
    #dt = 10800    # 3H -> sec 

    dy,dx = calc_dx_dy(lat[::-1],lon) # y,x grid spacing in meters
    dx = dx[np.newaxis,np.newaxis,:,:]
    dy = dy[np.newaxis,np.newaxis,:,np.newaxis]

    dP = np.gradient(lev)
    dP = dP[np.newaxis,:,np.newaxis,np.newaxis]

    ds = xr.open_dataset(f)
    
    T = ds['T'].values
    U = ds['U'].values
    V = ds['V'].values
    OMEGA = ds['OMEGA'].values
    theta = T*(100000/lev[np.newaxis,:,np.newaxis,np.newaxis])**(R/Cp)

    #dthetadt = np.gradient(theta,axis=0) / dt
    dthetadx = np.gradient(theta, axis=-1) / dx
    dthetady = np.gradient(theta,axis=-2)/ dy
    dthetadp = np.gradient(theta,axis=1)/ dP
    #Q = (T/theta)*(dthetadt + (U*dthetadx) + (V*dthetady) + (omega*dthetadp))
    Q = (T/theta)*((U*dthetadx) + (V*dthetady) + (OMEGA*dthetadp))
    Q = np.nan_to_num(Q,nan=0.0)
    
    Qkat = (T/theta)*((V*dthetady) + (OMEGA*dthetadp))
    Qkat = np.nan_to_num(Qkat,nan=0.0)
                   
    Q = Q.squeeze()
    theta = theta.squeeze()
    Qkat = Qkat.squeeze()
    
    # vertically integrated diabatic heating
    viQ = np.trapz(Q[:12],x=lev[:12,np.newaxis,np.newaxis],axis=0)
    
    return theta, Q, viQ, Qkat

def calc_divQ(Q, lat, lon, dx, dy):
    
    """
    Calculate divergence of diabatic heating rate.

    Parameters:
    Q (dtype: numpy.ndarray): 3D array with diabatic heating rate along pressure levels.
    lat (dtype: numpy.ndarray): 1D array of latitudes.
    lon (dtype: numpy.ndarray): 1D array of longitudes.
    dx (dtype: numpy.ndarray) : 2D array of grid spacing in the x-direction
    dy (dtype: numpy.ndarray) : 2D array of grid spacing in the y-direction

    Returns:
    divQ, dQdx, dQdy (dtype: tuple): divQ (divergence of Q), dQdx (x-component of divQ), 
                                    dQdy (y-component of divQ).
    """

    # Calculate gradients
    dQdx = (np.roll(Q, -1, axis=-1) - np.roll(Q, 1, axis=-1)) / (2 * dx)
    dQdy = (np.roll(Q, -1, axis=1) - np.roll(Q, 1, axis=1)) / (2 * dy)
    
    dQdx = dQdx.squeeze() * 3600
    dQdy = dQdy.squeeze() * 3600
    
    # Compute divergence
    divQ = dQdx + dQdy
    
    divQ = divQ.squeeze()

    return divQ,dQdx,dQdy

def read_cyclone_data(dates, hours, directory):
    
    """
    Read in cyclone track data as determined via use of CyTRACK from Alarcon et al., 2024 (https://github.com/apalarcon/CyTRACK).

    Parameters:
    dates (dtype: list): list of datetime objects as specified by arguments input file
    hours (dtype: list): list of 3H slices 
    directory (dtype: str): name of directory with surface data files.

    Returns:
    df (dtype: pandas.dataframe): pandas dataframe containing cyclone track data (YYYYMMDD, HH, MSLP, etc).
    """
    
    from datetime import datetime, timedelta
    
    # Combine the date and time from the first elements in dates and hours
    start_datetime_str = dates[0].strip() + hours[0].strip()
    start_datetime = datetime.strptime(start_datetime_str, '%Y%m%d%H')

    # Calculate the date 7 days before
    search_date = start_datetime# - timedelta(days=7)

    # Format the search date to match the file naming pattern
    search_date_str = search_date.strftime('%Y%m%d%H')

    # Search for the file in the given directory
    files = os.listdir(directory)
    target_file = None
    for file in files:
        if search_date_str in file:
            target_file = file
            break

    if not target_file:
        raise FileNotFoundError(f"No file found with date {search_date_str} in directory {directory}")

    file_path = os.path.join(directory, target_file)

    # Read the data from the found file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    data = []

    for line in lines:
        if line.strip().startswith("CySH"):
            # Skip lines starting with "CySH"
            continue

        parts = line.strip().split(',')
        if len(parts) > 2:  # To ensure it's a data line
            timestep = parts[0].strip() + '_' + parts[1].strip()  # Combining date and time to match timestep format
            date = parts[0]
            time = parts[1]
            lat = float(parts[2])
            lon = float(parts[3])
            slp = float(parts[4])
            data.append([timestep, date, time, lat, lon, slp])

    # Convert to DataFrame
    df = pd.DataFrame(data, columns=['timestep', 'date', 'time', 'lat', 'lon', 'slp'])
    return df

def plot_AR_detections(ar_data, vivt, fdate, latg, long, lon, lat, dates, hours, plot_cyclones): 
    
    """
    Plot polar projection maps of vIVT to indicate detected ARs in Southern Ocean; 
    additionally SLP field and detected cyclones are plotted correspondingly.

    Parameters:
    ar_data (dtype: list): lists of detected valid AR indices
    vivt (dtype: numpy.ndarray): 2D array of vertically integrated meridional vapor transport values 
    fdate (dtype: str): file datetime (YYYYMMDD-HH).
    latg (dtype: numpy.ndarray): 2D array of gridded latitude values 
    long (dtype: numpy.ndarray): 2D array of gridded longitude values
    lat (dtype: numpy.ndarray): 1D array of latitudes.
    lon (dtype: numpy.ndarray): 1D array of longitudes. 
    dates (dtype: list): list of datetime objects as specified by arguments input file
    hours (dtype: list): list of 3H slices 
    plot_cyclones (dtype: bool): True/False boolean that specifies 
                                whether cyclone data will be plotted concurrently

    Returns:
    landfall_lons (dtype: list): Where AR makes landfall with coastline 
                                according to mask specified by user input.
    """
    
    lon = lon - 180.
    
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import matplotlib.colors as mcolors
    
    if plot_cyclones: 
        cyc_dat = read_cyclone_data(dates,hours,'D:/Research/DATA/CyTRACK/CyTRACK_output')

    norm = plt.Normalize(0, 1500)
    cmap = mcolors.LinearSegmentedColormap.from_list("", ['white', 'lightsteelblue', 'cornflowerblue', 'royalblue', 'blue', 'limegreen', 'green', 'khaki', 'gold', 'orange', 'orangered', 'firebrick', 'darkorchid'])

    # Create a boolean mask initialized to True
    ar_mask = np.ones_like(vivt, dtype=bool)

    # Set the mask to False at locations corresponding to ar_data
    for j, i in ar_data:
        lat_idx = np.abs(latg[:, 0] - j).argmin()
        lon_idx = np.abs(long[0, :] - i).argmin()
        ar_mask[lat_idx, lon_idx] = False

    ar_mask = ar_mask.astype(int)
    
    coast_mask = np.array(nc.Dataset('D:/Research/DATA/masks/ANT_coast_mask.nc')['ANT_coast_mask'])
    coast_mask = coast_mask[::-1,:].astype(int)
    
    landfall_mask = np.logical_and(ar_mask, coast_mask)

    # Step 2: Extract Longitudes where Overlap Occurs
    landfall_indices = np.where(landfall_mask)
    landfall_lons = np.unique(landfall_indices[1])  # Assuming longitude is the second axis in the mask
    
    ds_slp = xr.open_dataset('D:/Research/DATA/MERRA2/M2T1NXSLV.5.12.4/3H_slices/{}.nc'.format(fdate))
    slp = ds_slp['SLP'].values / 100.

    fig = plt.figure(figsize=(12, 12))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180, 180, -90, -40], ccrs.PlateCarree())

    ax.coastlines(resolution="110m", color='k', linewidth=3, zorder=3)
    gl = ax.gridlines(linestyle='--', color='grey', draw_labels=True, alpha=0.7,zorder=2)
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 10))
    gl.ylocator = plt.FixedLocator(np.arange(-90, -40, 10))

    lvls = np.linspace(0, 1500, 21)
    cf = plt.contourf(lon, lat, vivt, cmap=cmap, norm=norm, levels=lvls, extend='both', transform=ccrs.PlateCarree(), alpha=0.9)

    #c = plt.contour(lon, lat, ar_mask, levels=[0.5], colors=['k'], linewidths=2, linestyles='--', transform=ccrs.PlateCarree(), zorder=10)
    c = plt.contour(lon, lat, slp[::-1,:], levels=np.arange(940,1040,10), colors=['k'], linewidths=2, linestyles='--', transform=ccrs.PlateCarree(), zorder=10)
    plt.clabel(c, inline=True, fontsize=14)

    cb = plt.colorbar(cf, orientation="vertical", extend='both', fraction=0.046, pad=0.02)
    cb.set_label('vIVT (kg m-1 s-1)', size=8, rotation=90, labelpad=12)
    cb.ax.tick_params(labelsize=10)
    
    # Plot cyclone centers and pressures for the current fdate
    if plot_cyclones:
        cycs = cyc_dat[cyc_dat['timestep'] == fdate]
        for idx, row in cycs.iterrows():
            ax.scatter(row['lon'], row['lat'], color='red', s=65, marker='*', transform=ccrs.PlateCarree(), zorder=5)
            ax.text(row['lon'], row['lat'], f"{row['slp']:.1f}", color='firebrick', fontsize=14, ha='right', transform=ccrs.PlateCarree(), zorder=6)
    
    plt.title('AR detections - {}'.format(fdate))
    plt.savefig('D:/Research/DATA/wille_ar_detections/plots/ARs_cyclones/AR_detections_{}.png'.format(fdate))
    plt.close()
    
    return landfall_lons

def plot_diabatic_heating(viQ,Q,Qkat,theta,vivt,U,V,dx,dy,fdate,lat,lon,lev,dates,hours,landfall_lons,plot_cyclones):
    
    """
    Plot polar projection maps of diabatic heating due to cyclones/ARs detected ARs in Southern Ocean.

    Parameters:
    viQ (dtype: numpy.ndarray): 2D array of vertically integrated diabatic heating (1000hPa - 700hPa (surface - ~2500m)).
    Q (dtype: numpy.ndarray): 3D array of diabatic heating rate.
    Qkat (dtype: str): 3D array of x-component of diabatic heating rate.
    theta (dtype: numpy.ndarray): 3D array of potential temperature. 
    vivt (dtype: numpy.ndarray):  2D array of vertically integrated meridional vapor transport values (1000hPa - 200hPa).
    U (dtype: numpy.ndarray): 2D array of eastward wind vectors.
    V (dtype: numpy.ndarray): 2D array of northward wind vectors.
    dx (dtype: numpy.ndarray) : 2D array of grid spacing in the x-direction.
    dy (dtype: numpy.ndarray) : 2D array of grid spacing in the y-direction.
    fdate (dtype: str): file datetime (YYYYMMDD-HH).
    lat (dtype: numpy.ndarray): 1D array of latitudes.
    lon (dtype: numpy.ndarray): 1D array of longitudes.
    lev (dtype: numpy.ndarray) : 1D array of pressure levels in Pa.
    dates (dtype: list): list of datetime objects as specified by arguments input file.
    hours (dtype: list): list of 3H slices. 
    landfall_lons (dtype: list): Where AR makes landfall with coastline 
                                according to mask specified by user input.
    plot_cyclones (dtype: bool): True/False boolean that specifies 
                                whether cyclone data will be plotted concurrently.

    Returns:
    None (dtype: None): Saves .png files out to directory specified by user in this function.
    """
    
    lon = lon - 180.
    
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    
    if plot_cyclones: 
        cyc_dat = read_cyclone_data(dates,hours,'D:/Research/DATA/CyTRACK/CyTRACK_output') 
        
    ds_surf = xr.open_dataset('D:/Research/DATA/MERRA2/M2T1NXSLV.5.12.4/3H_slices/{}.nc'.format(fdate))
    slp = ds_surf['SLP'].values / 100.
    v_surf = ds_surf['V10M'].values
    u_surf = ds_surf['U10M'].values
    
    vQkat = np.trapz(Qkat[:12],x=lev[:12,np.newaxis,np.newaxis],axis=0)
    
    ######### PLOT #1 : Diabatic Heating along the coastline of the Antarctic continent #########
    
    fig = plt.figure(figsize=(12, 12))
    ax = plt.axes(projection=ccrs.SouthPolarStereo())
    ax.set_extent([-180, 180, -90, -55], ccrs.PlateCarree())

    ax.coastlines(resolution="110m", color='k', linewidth=3, zorder=4)
    gl = ax.gridlines(linestyle='--', color='w', draw_labels=True)
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 10))
    gl.ylocator = plt.FixedLocator(np.arange(-90, -54, 10))

    ax.quiver(lon[::5], lat[::5], v_surf[::-1,:][::5,::5], v_surf[::-1,:][::5,::5], transform=ccrs.PlateCarree(), color='k', width= 0.0015, zorder=6)
    
    c = plt.contour(lon, lat, slp[::-1,:], levels=np.arange(940,1040,10), colors=['k'], linewidths=2, linestyles='--', transform=ccrs.PlateCarree(), zorder=10)
    plt.clabel(c, inline=True, fontsize=14)

    #lvls = np.arange(-5,5,0.5)
    #cf = plt.contourf(lon, lat, vQkat[::-1,:], cmap=plt.cm.seismic_r, levels=lvls, extend = 'both', transform=ccrs.PlateCarree(), zorder=2)
    lvls = np.arange(-10,10+1,1)
    cf = plt.contourf(lon, lat, viQ[::-1,:], cmap=plt.cm.seismic, levels=lvls, extend = 'both', transform=ccrs.PlateCarree(), zorder=2)

    # Plot cyclone centers and pressures for the current fdate
    if plot_cyclones:
        cycs = cyc_dat[cyc_dat['timestep'] == fdate]
        for idx, row in cycs.iterrows():
            ax.scatter(row['lon'], row['lat'], color='k', s=90, marker='*', transform=ccrs.PlateCarree(), zorder=5)
            ax.text(row['lon'], row['lat'], f"{row['slp']:.1f}", color='k', fontsize=14, ha='right', weight='bold', transform=ccrs.PlateCarree(), zorder=6)

    cb = plt.colorbar(cf, orientation="vertical", extend='both', fraction=0.046, pad=0.02)
    #cb.set_label('Diabatic Heating Rate (K hr-1)',size=8,rotation=90,labelpad=12)
    cb.set_label('Diabatic Heating Rate (K kg m-2 s-1)',size=8,rotation=90,labelpad=12)
    cb.ax.tick_params(labelsize=10)
    #plt.title('Q @ 700hPa (K hr-1) {}UTC'.format(fdate))
    plt.title('vertically integrated Q (K kg m-2 s-1) {}UTC'.format(fdate))
    #plt.show()
    plt.savefig('D:/Research/DATA/wille_ar_detections/plots/Q_AR_cyclones/Q_AR_detections_{}.png'.format(fdate))
    plt.close()
    
    ######### PLOT #2 : Diabatic Heating along columns of the atmosphere near landfalling AR #########
    '''
    import matplotlib.colors as mcolors
    
    #dP = np.gradient(lev)
    
    # Step 3: Find contiguous groups of longitudes
    contiguous_groups = find_contiguous_groups(landfall_lons)

    if len(contiguous_groups) > 0:
        # Define levels
        levels = np.linspace(-2.5, 2.5, 21)  # Change levels as needed

        # Create a normalizer
        norm = mcolors.BoundaryNorm(boundaries=levels, ncolors=256)

        # Create a colormap and set masked values to be black
        cmap = plt.cm.viridis
        cmap.set_bad(color='grey')

        for group in contiguous_groups:
            start = np.abs(lon - group[0]).argmin()
            stop = np.abs(lon - group[-1]).argmin() + 1  # Inclusive of stop

            # Ensure start and stop are within bounds
            start = max(0, start)
            stop = min(len(lon), stop)

            if start >= stop:
                continue

            plt.figure(figsize=(14, 8))

            # Mask invalid values in Q
            Q_mean_masked = np.ma.masked_invalid(np.mean(Q[::-1, :, start:stop], axis=-1) * 3600)

            pcm = plt.pcolormesh(lat, lev[::-1] / 100, Q_mean_masked, cmap=cmap, norm=norm)

            quiv = plt.quiver(lat, lev[::-1] / 100, np.mean(U[::-1, :, start:stop], axis=-1), 
                              np.mean(V[::-1, :, start:stop], axis=-1), color='k', width=0.001)

            cont = plt.contour(lat, lev[::-1] / 100, np.mean(theta[::-1, :, start:stop], axis=-1),
                               levels=np.linspace(250, 350, 21), colors='k', linestyles='--', linewidths=1)
            
            
            
            plt.clabel(cont, inline=True, fontsize=12)
            plt.gca().invert_yaxis()

            # Plot cyclone centers
            cycs_filtered = cycs[(cycs['lon'] >= group[0]) & (cycs['lon'] <= group[-1])]
            for idx, row in cycs_filtered.iterrows():
                plt.scatter(row['lat'], lev[8] / 100, color='limegreen', s=90, marker='*')
                plt.text(row['lat'], lev[8] / 100, f"{row['slp']:.1f}", color='limegreen', fontsize=14, ha='right', weight='bold')

            cb = plt.colorbar(pcm)
            cb.set_label('Q (K hr-1)', size=8, rotation=90, labelpad=12)
            plt.title('Diabatic Heating Contribution {}UTC (lon {} to {})'.format(fdate, lon[start], lon[stop-1]))
            plt.savefig('D:/Research/DATA/wille_ar_detections/plots/Q_landfall/Q_landfall_{}_LON_{}-{}.png'.format(fdate,lon[start],lon[stop-1]))
            plt.close()
    else: 
        pass
    '''
    
def plot_vdQdy_cycs_ARs(Q,vivt,imfc,lat,lon,dx,dy,lev,ar_data,plot_cyclones,fdate,dates,hours):
    
    """
    Plot polar projection maps of vertically integrated diabatic heating 
    (MERIDIONAL COMPONENT ONLY) due to cyclones/ARs detected ARs in Southern Ocean.

    Parameters:
    Q (dtype: numpy.ndarray): 3D array of diabatic heating rate.
    vivt (dtype: numpy.ndarray):  2D array of vertically integrated meridional vapor transport values (1000hPa - 200hPa).
    imfc (dtype: numpy.ndarray):  2D array of vertically integrated moisture flux convergence (1000hPa - 700hPa).
    lat (dtype: numpy.ndarray): 1D array of latitudes.
    lon (dtype: numpy.ndarray): 1D array of longitudes.
    dx (dtype: numpy.ndarray) : 2D array of grid spacing in the x-direction.
    dy (dtype: numpy.ndarray) : 2D array of grid spacing in the y-direction.
    lev (dtype: numpy.ndarray) : 1D array of pressure levels in Pa.
    ar_data (dtype: list): lists of detected valid AR indices
    plot_cyclones (dtype: bool): True/False boolean that specifies 
                                whether cyclone data will be plotted concurrently.
    fdate (dtype: str): file datetime (YYYYMMDD-HH).
    dates (dtype: list): list of datetime objects as specified by arguments input file.
    hours (dtype: list): list of 3H slices. 

    Returns:
    None (dtype: None): Saves .png files out to directory specified by user in this function.
    """
    
    import matplotlib.pyplot as plt
    import cartopy.crs as ccrs
    import matplotlib.colors as mcolors
    
    norm = plt.Normalize(-0.5,0.5+0.1)
    cmap = mcolors.LinearSegmentedColormap.from_list("", ['lightsteelblue', 'cornflowerblue', 'royalblue', 'blue', 'green', 'limegreen', 'white', 'khaki', 'gold', 'orange', 'orangered', 'firebrick', 'darkorchid'])

    lat = lat[::-1]
    lon = lon - 180.
    
    # Create a boolean mask initialized to True
    ar_mask = np.ones_like(vivt, dtype=bool)

    # Set the mask to False at locations corresponding to ar_data
    for j, i in ar_data:
        lat_idx = np.abs(lat - j).argmin()
        lon_idx = np.abs(lon - i).argmin()
        ar_mask[lat_idx, lon_idx] = False

    ar_mask = ar_mask.astype(int)
    
    if plot_cyclones: 
        cyc_dat = read_cyclone_data(dates,hours,'D:/Research/DATA/CyTRACK/CyTRACK_output')
    
    divQ,dQdx,dQdy = calc_divQ(Q,lat,lon,dx,dy)

    vdQdy = np.trapz(dQdy[:17],x=lev[:17,np.newaxis,np.newaxis],axis=0)

    #levels=np.arange(-0.5,0.5+0.05,0.05)
    levels=np.arange(-50,50+10,10)

    fig = plt.figure(figsize=(20, 20))

    # Define the projection
    projection = ccrs.SouthPolarStereo()

    # Create subplots
    ax1 = plt.axes(projection=ccrs.SouthPolarStereo())

    # Set extent
    ax1.set_extent([-180, 180, -90, -50], ccrs.PlateCarree())

    ax1.coastlines(resolution="110m", color='k', linewidth=3, zorder=3)
    gl = ax1.gridlines(linestyle='--', color='grey', draw_labels=True, alpha=0.7, zorder=2)
    gl.xlocator = plt.FixedLocator(np.arange(-180, 181, 10))
    gl.ylocator = plt.FixedLocator(np.arange(-90, -40, 10))
    gl.xlabel_style = {'size': 10}
    gl.ylabel_style = {'size': 10}
    
    #c = plt.contour(lon, lat, ar_mask, levels=[0.5], colors=['k'], linewidths=2, linestyles='--', transform=ccrs.PlateCarree(), zorder=10)
    # Plot the meridional gradient in diabatic heating 
    c1 = ax1.contourf(lon, lat, vdQdy, cmap=cmap, levels=levels, extend = 'both', transform=ccrs.PlateCarree())
    #c1 = ax1.contourf(lon, lat, imfc*3600, cmap=cmap, levels=levels, extend = 'both', transform=ccrs.PlateCarree())
    if plot_cyclones:
        cycs = cyc_dat[cyc_dat['timestep'] == fdate]
        for idx, row in cycs.iterrows():
            ax1.scatter(row['lon'], row['lat'], color='k', s=90, marker='*', transform=ccrs.PlateCarree(), zorder=5)
            ax1.text(row['lon'], row['lat'], f"{row['slp']:.1f}", color='k', fontsize=14, ha='right', weight='bold', transform=ccrs.PlateCarree(), zorder=6)

    cb = plt.colorbar(c1, orientation="vertical", extend='both', fraction=0.046, pad=0.02)
    #cb.set_label('Diabatic Heating Rate (K hr-1)',size=8,rotation=90,labelpad=12)
    cb.set_label('K Pa m-1 hr-1',size=8,rotation=90,labelpad=12)
    cb.ax.tick_params(labelsize=10)
    #plt.title('Q @ 700hPa (K hr-1) {}UTC'.format(fdate))
    plt.title('Vertically Integrated Meridional Gradient of Q (K Pa m-1 s-1) {}UTC'.format(fdate))
    #plt.show()
    plt.savefig('D:/Research/DATA/wille_ar_detections/plots/vdQdy_cyclones_ARs/vdQdy_cyclones_ARs_{}.png'.format(fdate))
    #plt.savefig('D:/Research/DATA/wille_ar_detections/plots/imfc_cyclones_ARs/imfc_cyclones_ARs_{}.png'.format(fdate))
    plt.close()
    
def wille_vivt_ar_detect(idir,dates,hours,plot_cyclones, plot_Q):
    
    """
    Runs vIVT AR detection scheme from Wille et al., 2019 (https://github.com/jonathanwille/Antarctic-lab/tree/v2.3).

    Parameters:
    idir (dtype: list): directory where data files are located
    dates (dtype: list): list of datetime objects as specified by arguments input file.
    hours (dtype: list): list of 3H slices. 
    plot_cyclones (dtype: bool): True/False boolean that specifies 
                                whether cyclone data will be plotted concurrently.
    plot_Q (dtype: bool): True/False boolean that specifies 
                                whether diabatic heating data will be plotted concurrently.

    Returns:
    None (dtype: None): Saves out .txt files with detected AR lat/lons, as well as 
    .png files using plotting functions defined above.
    """

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
    lev = ex_ds.variables['lev'][:] * 100    # hPa --> Pa
    
    lat = lat[::-1]
    lon = lon + 180

    long, latg = np.meshgrid(lon, lat)
    
    coast_mask = np.array(nc.Dataset('D:/Research/DATA/masks/ANT_coast_mask.nc')['ANT_coast_mask'])
    coast_mask = coast_mask[::-1,:].astype(bool)
    
    g0 = 9.80665  # Gravitational acceleration in m/s^2
    #vivt = np.zeros((len(flist)*nt,nlat,nlon))
    
    dy,dx = calc_dx_dy(lat[::-1],lon) # y,x grid spacing in meters
    dx = dx[np.newaxis,np.newaxis,:,:]
    dy = dy[np.newaxis,np.newaxis,:,np.newaxis]

    # Iterate over daily files
    for f in flist:
        
        fdate = os.path.basename(f).split('.')[0]
        
        print("   ---> | Processing : " + f)
        
        prev_mon_files = get_previous_30_days_files(idir,dates,fdate)
        thresh = calculate_threshold_vivt(prev_mon_files, lat, lon, lev,g0)
        
        lat_index_list = []
        lon_index_list = []
    #for f in flist[:1]:
        # read the file in and extract daily 3H vIVT
        with xr.open_dataset(f) as ds:
            
            v = ds['V'].values
            u = ds['U'].values
            qv = ds['QV'].values

            qvv = np.trapz(qv * v, x=lev[:, np.newaxis, np.newaxis], axis=0)
            #qvu = np.trapz(qv * u, x=lev[:, np.newaxis, np.newaxis], axis=0)
            #vivt = (1 / g0) * np.sqrt(qvu**2+qvv*2)
            vivt = (1/g0) * qvv

            vivt = vivt[::-1, :]
            vivt = np.ma.array(data=vivt,mask=coast_mask)
            
            # Calculate gradients
            duqdx = ((np.roll(u*qv, -1, axis=-1) - np.roll(u*qv, 1, axis=-1)) / (2 * dx)).squeeze()
            dvqdy = ((np.roll(v*qv, -1, axis=1) - np.roll(v*qv, 1, axis=1)) / (2 * dy)).squeeze()
            
            #duqdx = (np.gradient(u*qv,axis=-1)/dx).squeeze()
            #dvqdy = (np.gradient(v*qv,axis=1)/dy).squeeze()
            
            imfc = np.trapz(duqdx[:12]+dvqdy[:12], x=lev[:12, np.newaxis, np.newaxis], axis=0)
            
            if plot_Q:
                theta, Q, viQ, Qkat = calculate_diabatic_heating(f,lat,lon,lev)

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
            #ar_lon=lon-180
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
            
            # Create a boolean mask initialized to True
            ar_mask = np.ones_like(vivt, dtype=bool)

            # Set the mask to False at locations corresponding to ar_data
            for j, i in ar_data:
                lat_idx = np.abs(latg[:, 0] - j).argmin()
                lon_idx = np.abs(long[0, :] - i).argmin()
                ar_mask[lat_idx, lon_idx] = False

            ar_mask = ar_mask.astype(int)
            
            ar_landfall_pts = find_ar_landfall_points(ar_mask, lat, lon)
            landfall_lons = [i for _, i in ar_landfall_pts]
            
            plot_AR_detections(ar_data, vivt, fdate, latg, long, lon, lat, dates, hours, plot_cyclones)
            
            if plot_Q: 
                plot_diabatic_heating(viQ,Q,Qkat,theta,vivt,u[:,::-1,:],v[:,::-1,:],dx,dy,fdate,lat,lon,lev,dates,hours,landfall_lons,plot_cyclones)
                
            plot_vdQdy_cycs_ARs(Q,vivt,imfc,lat,lon,dx,dy,lev,ar_data,plot_cyclones,fdate,dates,hours)