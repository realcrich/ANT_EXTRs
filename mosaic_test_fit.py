# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 17:15:40 2024

@author: Collin
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from funcs.get_shp_coords import get_shp_coords

import h5py

ddir = 'D:/Research/DATA/test_fit'

class mosaic:
    """Utility for creating spatial mosaics
    """
    def __init__(self, **kwargs):
        self.extent = [np.inf,-np.inf,np.inf,-np.inf]
        self.spacing = [None,None]
        self.fill_value = np.nan

    def update_spacing(self, x, y):
        """
        update the step size of mosaic
        """
        try:
            self.spacing = (x[1] - x[0], y[1] - y[0])
        except:
            pass
        return self

    def update_bounds(self, x, y):
        """
        update the bounds of mosaic
        """
        # check that there is data
        if not np.any(x) or not np.any(y):
            return self
        # get extent of new data
        extent = [x.min(), x.max(), y.min(), y.max()]
        if (extent[0] < self.extent[0]):
            self.extent[0] = np.copy(extent[0])
        if (extent[1] > self.extent[1]):
            self.extent[1] = np.copy(extent[1])
        if (extent[2] < self.extent[2]):
            self.extent[2] = np.copy(extent[2])
        if (extent[3] > self.extent[3]):
            self.extent[3] = np.copy(extent[3])
        return self

    def image_coordinates(self, x, y):
        """
        get the image coordinates
        """
        # check that there is data
        if not np.any(x) or not np.any(y):
            return (None, None)
        # get the image coordinates
        iy = np.array((y[:,None] - self.extent[2])/self.spacing[1], dtype=np.int64)
        ix = np.array((x[None,:] - self.extent[0])/self.spacing[0], dtype=np.int64)
        return (iy, ix)

    @property
    def dimensions(self):
        """Dimensions of the mosaic"""
        dims = [None, None]
        # calculate y dimensions with new extents
        dims[0] = np.int64((self.extent[3] - self.extent[2])/self.spacing[1]) + 1
        # calculate x dimensions with new extents
        dims[1] = np.int64((self.extent[1] - self.extent[0])/self.spacing[0]) + 1
        return dims

    @property
    def shape(self):
        """Shape of the mosaic"""
        return (self.dimensions[0], self.dimensions[1], )

    @property
    def x(self):
        """X-coordinates of the mosaic"""
        return self.extent[0] + self.spacing[0]*np.arange(self.dimensions[1])

    @property
    def y(self):
        """Y-coordinates of the mosaic"""
        return self.extent[2] + self.spacing[1]*np.arange(self.dimensions[0])

# PURPOSE: read a variable group from ICESat-2 ATL15
def read_test_fit(infile, group='fit_statistics', fields=None):
    # dictionary with ATL15 variables
    d = {}
    
    with h5py.File(infile, 'r') as h5f:
        if fields is None:
            fields = ['h_fit', 'x', 'y', 'time']
        
        for field in fields:
            if field in h5f[group]:
                d[field] = h5f[group][field][:]
            else:
                print(f"Warning: Field '{field}' not found in group '{group}'")
    
    return d


flist=[]
for fp in os.listdir(ddir): 
    flist.append(os.path.join(ddir,fp))

def mosaic_test_fit(files, mosaic): 
    # create mosaic of ATL15 data
    # iterate over each ATL15 file
    for f in files:
        # get ATL15 dimension variables from group
        d = read_test_fit(f, group='fit_statistics', fields=['x','y','time'])
        # update the mosaic grid spacing
        mosaic.update_spacing(d['x'], d['y'])
        mosaic.update_bounds(d['x'], d['y'])

    # dimensions of output mosaic
    ny, nx = mosaic.shape
    nt = len(d['time'])
    print(f"Mosaic shape: ny={ny}, nx={nx}, nt={nt}")

    # create output mosaic
    h_fit = {}
    h_fit['x'] = np.copy(mosaic.x)
    h_fit['y'] = np.copy(mosaic.y)
    h_fit['time'] = np.copy(d['time'])
    valid_mask = np.zeros((nt,ny,nx), dtype=bool)

    # Initialize h_fit with the correct shape
    h_fit['h_fit'] = np.ma.zeros((nt,ny,nx), fill_value=np.nan)
    h_fit['h_fit'].mask = np.ones((nt,ny,nx), dtype=bool)

    # iterate over each ATL15 file
    for f in files:
        # get ATL15 variables from group
        d = read_test_fit(f, group='fit_statistics')
        # get the image coordinates of the input file
        iy, ix = mosaic.image_coordinates(d['x'], d['y'])
        print(f"Shape of iy: {iy.shape}, Shape of ix: {ix.shape}")
        print(f"Shape of d['h_fit']: {d['h_fit'].shape}")

        # Ensure iy and ix are within bounds
        valid = (iy >= 0) & (iy < ny) & (ix >= 0) & (ix < nx)
        print(f"Shape of valid: {valid.shape}")

        # Create a meshgrid of indices
        iy_mesh, ix_mesh = np.meshgrid(iy.flatten(), ix.flatten(), indexing='ij')

        # Get the number of time steps in this file
        nt_file = min(d['h_fit'].shape[2], nt)

        # Use numpy's advanced indexing
        valid_mask[:nt_file, iy_mesh, ix_mesh] |= valid

        # Assign values
        h_fit_file = d['h_fit'].transpose(2,0,1)[:nt_file]
        h_fit['h_fit'][:nt_file, iy_mesh, ix_mesh] = np.where(
            valid[np.newaxis, :, :],
            h_fit_file,
            h_fit['h_fit'][:nt_file, iy_mesh, ix_mesh]
        )

    # update masks for variables
    h_fit['h_fit'].mask = np.logical_not(valid_mask)

    return h_fit

mosaic=mosaic()
# Call the function
h_fit = mosaic_test_fit(flist, mosaic)

x_grid, y_grid = np.meshgrid(h_fit['x'],h_fit['y'])

polygons = get_shp_coords('ANT')

for yr in range(0, h_fit['h_fit'].shape[0]):
    plt.figure(figsize=(10,8))

    for polygon in polygons:
        xi,yi = polygon
        plt.plot(xi, yi,color='k',lw=1.5, zorder=4)
    #plt.plot(xint, yint, 'm', lw = 2.5)
    levels = np.linspace(-3000, 3000, 21)
    res = plt.contourf(x_grid, y_grid, h_fit['h_fit'][yr, :, :], levels = levels, cmap='RdBu', extend='both',zorder=1)

# Add contour lines
    #contours = plt.contour(x_grid, y_grid, dplt[yr, :, :], levels=levels, colors='k', linewidths=0.35)
    #plt.clabel(contours, inline=True, fontsize=7)

    #plt.xlim(x_min, x_max)
    #plt.ylim(y_min, y_max)

    #plt.xlim(-.9e6, 2.35e6)
    #plt.ylim(.75e6, 2.25e6)

    #x_min, x_max = np.min(x_grid), np.max(x_grid)
    #y_min, y_max = np.min(y_grid), np.max(y_grid)

    #plt.xlim(-0.9e6, x_max)
    #plt.ylim(y_min+0.5e6, y_max-0.5e6)

    plt.grid()
    plt.colorbar(res, orientation='vertical', label='dhdt', extend='both', fraction=0.046, pad=0.04)

    plt.gca().set_xticks([])  # Remove x-axis tick labels
    plt.gca().set_yticks([])  # Remove y-axis tick labels

    plt.title('{}'.format(2002+yr))