# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 11:06:33 2024

@author: Collin
"""

# make_poly reads in a shapefile_path for Antarctic IMBIE Rignot Drainage Basins
# creating a shapely Polygon from the .shp file and (x,y) coordinate arrays for the exterior of the polygon
# also returns name of the selected region 

# reg_choice = user basin of interest to create polygon for (integer value)
# interior = additional 'y/n' option to read in and create polygon for East Antarctic interior (above 2500m elevation)

from PIL import Image
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
import shapefile
import numpy as np
from funcs.map_ll import map_ll

def make_poly(reg_choice, interior):
    # Prompt the user for the path to the PNG image
    png_path = 'D:/Research/drainage_basins.jpg'

    # Display the PNG image inline using matplotlib
    regions_map = Image.open(png_path)
    plt.imshow(regions_map)
    plt.axis('off')
    plt.show()

    # Prompt the user for the path to the shapefile
    shapefile_path = 'D:/Research/DATA/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'

    shape = shapefile.Reader(shapefile_path)
    records = shape.records()
    '''
    # Print the choices for the reg variable based on the shapefile records
    print("Available choices for 'reg':")
    for i, record in enumerate(records):
        print(f"{i}: {record}")  # Adjust this line based on the actual attribute field you want to display
    '''
    # Create a list of basin names
    basin_names = [
        'islands', 'hhp', 'fg', 'eep', 'ddp', 'cpd', 'bc', 'aap', 'jppk', 'gh',
        'dpe', 'apb', 'ccp', 'ka', 'jjpp', 'ippj', 'iipp', 'hpi', 'epf'
    ]
    
    # Check if reg_choice is a valid basin name
    if reg_choice not in basin_names:
        raise ValueError(f"Invalid basin name: {reg_choice}")

    # Get the index of the selected basin name
    basin_index = basin_names.index(reg_choice)

    rec = shape.shapeRecords()[basin_index]
    points = rec.shape.points

    x = np.array(points)[:, 0]
    y = np.array(points)[:, 1]

    reg_coords = np.array(list(zip(x, y)))

    # Create the polygon using the generated name
    poly_reg = Polygon(reg_coords)
    poly_reg.name = reg_choice
    
    if interior == 'y':
        # Prompt the user for the file path
        ant2500_shp_path = 'D:/Research/DATA/eais_2500m_elevation.ascii'

        # Open the file
        f = open(ant2500_shp_path, 'r')
        data = np.genfromtxt(f)
        x2500, y2500 = map_ll(data[:, 0], data[:, 1], HEM='S', FORMAT='tuple')
        grid_pts_2500 = np.array(list(zip(x2500, y2500)))
        poly_2500 = Polygon(grid_pts_2500)
        
        poly_coast = poly_reg.difference(poly_2500)
    
        return x, y, poly_reg.name, poly_reg, x2500, y2500, poly_2500, poly_coast
    
    else: 
        
        return x, y, poly_reg.name, poly_reg, None, None, None, None