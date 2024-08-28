# -*- coding: utf-8 -*-
"""
Created on Mon Jan  8 14:04:24 2024

@author: Collin
"""

# reads in polygon and associated coordinates and makes limited x and y arrays 
# based on the region for plotting extents 

def limit_data_bbox(x, y, poly, x_length, y_length):
    
    x_poly, y_poly = poly.exterior.xy
    
    len_x = max(x_poly)-min(x_poly) + 0.3e6
    len_y = max(y_poly)-min(y_poly) + 0.3e6
    
    # Calculate the center point of the polygon
    center_x = (min(x_poly) + max(x_poly)) / 2
    center_y = (min(y_poly) + max(y_poly)) / 2
    
    # Calculate the coordinates of the bounding box corners
    x_min = center_x - len_x / 2
    x_max = center_x + len_x / 2
    y_min = center_y - len_y / 2
    y_max = center_y + len_y / 2

    return x_min , x_max, y_min, y_max