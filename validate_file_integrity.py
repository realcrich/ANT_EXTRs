# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 11:09:24 2024

@author: Collin
"""

import os
import xarray as xr

def check_file_integrity(file_path):
    """
    Check the integrity of a NetCDF file by attempting to open it and read a variable.
    """
    try:
        with xr.open_dataset(file_path, engine="netcdf4") as ds:
            # Attempt to read a small subset of data from each variable
            for var in ds.data_vars:
                _ = ds[var].values[0, 0, 0]  # Read a single value
        return True, None
    except Exception as e:
        return False, str(e)

def check_directory_integrity(directory_path):
    """
    Check the integrity of all NetCDF files in a directory.
    """
    corrupted_files = []
    all_files = [os.path.join(directory_path, fp) for fp in os.listdir(directory_path) if fp.endswith('.nc')]
    all_files = sorted(all_files)

    for file_path in all_files:
        is_valid, error_message = check_file_integrity(file_path)
        if not is_valid:
            corrupted_files.append((file_path, error_message))

    return corrupted_files

# Define the directory path
dp = 'D:/Research/DATA/MERRA2/M2I3NPASM.5.12.4/3H_slices'

# Check the integrity of all files in the directory
corrupted_files = check_directory_integrity(dp)

# Report the results
if corrupted_files:
    print("The following files are corrupted:")
    for file_path, error in corrupted_files:
        print(f"File: {file_path} - Error: {error}")
else:
    print("All files are intact.")