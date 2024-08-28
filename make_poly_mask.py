# -*- coding: utf-8 -*-
"""
Created on Sun Jan  7 18:23:59 2024

@author: Collin
"""

import numpy as np
import shapely.geometry as geometry

def make_poly_mask(data, poly, x_grid, y_grid):
    x_poly, y_poly = poly.exterior.xy

    # Create a mask array with True values
    mask = np.ones(data.shape, dtype=bool)

    # Create a MultiPoint object from the coordinates
    points = geometry.MultiPoint(np.column_stack((x_grid.flatten(), y_grid.flatten())))

    # Check if the points are within the polygon
    within_polygon = np.array([poly.contains(point) for point in points])

    # Stack the within_polygon array to repeat it as many times as the length of data
    within_polygon_stacked = np.stack([within_polygon] * data.shape[0], axis=0)

    # Reshape within_polygon_stacked to match the shape of mask
    within_polygon_stacked = within_polygon_stacked.reshape(mask.shape)

    # Update the mask based on the within_polygon values
    mask[within_polygon_stacked] = False

    return x_poly, y_poly, mask