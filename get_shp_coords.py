# -*- coding: utf-8 -*-
"""
Created on Fri Jan  5 13:41:10 2024

@author: Collin
"""

# read in .shp file and extract all polygon geoms to create list of all points for plotting

import geopandas as gpd

def get_shp_coords(reg):
    # Prompt the user for the path to the shapefile
    if reg == 'ANT':
        shp_path = 'D:/Research/DATA/ANT_Basins_IMBIE2_v1.6/ANT_Basins_IMBIE2_v1.6.shp'
    else:
        pass
    # Read the shapefile using geopandas
    gdf = gpd.read_file(shp_path)

    # Extract x, y coordinates for each polygon
    polygons = []
    for idx, row in gdf.iterrows():
        geometry = row['geometry']
        if geometry.geom_type == 'Polygon':
            exterior_coords = list(geometry.exterior.coords)
            x, y = zip(*exterior_coords)
            polygons.append((x, y))
        elif geometry.geom_type == 'MultiPolygon':
            for polygon in geometry:
                exterior_coords = list(polygon.exterior.coords)
                x, y = zip(*exterior_coords)
                polygons.append((x, y))

    return polygons