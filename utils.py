"""
Created by:
Riddhima Puri

Utility functions used by urmask.

This module contains functions adapted from the URCLIMASK package developed by Diez-Sierra et al. (2025).

Original source:
Diez-Sierra, J., Quintana, Y., Langendijk, G.S. et al. (2025).
A global CORDEX-based dataset delineating urban areas and their surroundings to assess climate change in megacities.
Scientific Data, 12, 1961.
https://doi.org/10.1038/s41597-025-06257-1

Modifications have been made where necessary to support the UrbanRuralSelection framework and CPRCM applications.
"""

import numpy as np
import xarray as xr
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon
from shapely.ops import unary_union

def find_city_coords(ds: xr.DataArray, lat_city: float, lon_city: float):
    """
    Returns closest grid cell to a city point e.g. [lat_city, lon_city] in geographic coordinates.

    Parameters:
    -----------
    ds : xarray.Dataset
        The dataset containing latitude and longitude coordinates.

    Returns:
    ---------
    lat_city_ds, lon_city_ds : float
        Latitude and longitude of the closest grid cell to the city point.
    ilat, ilon : int
        Indices of the closest grid cell to the city point in the dataset.
    """
    # Calculate distance from each grid cell to the city point
    dist = (ds['lon']-lon_city)**2 + (ds['lat']-lat_city)**2
    # Get indices of the closest grid cell to the city point
    ilat, ilon = np.unravel_index(dist.argmin(), dist.shape)
    # Latitude and longitude of the closest grid cell to the city point
    lat_city_ds = ds['lat'][ilat, ilon].values
    lon_city_ds = ds['lon'][ilat, ilon].values
    
    return lat_city_ds, lon_city_ds, ilat, ilon


def plot_urban_polygon(mask, ax):
    '''
    Plots urban and non-urban polygons from a mask dataset on the given axis.
    Returns GeoDataFrames for urban and non-urban areas.
    '''
    # Assume the mask is in the 'urmask' variable
    lon2d = mask.lon.values
    lat2d = mask.lat.values
    
    # Create lists to store polygons for urban areas (mask == 1) and non-urban areas (mask == 0)
    urban_polygons = []
    non_urban_polygons = []
    
    if lon2d.ndim == 1:
        dist_lat = abs(lat2d[1] - lat2d[0]) / 2
        dist_lon = abs(lon2d[1] - lon2d[0]) / 2
        for lon in range(len(lon2d)):
            for lat in range(len(lat2d)):
                square = Polygon([
                    (round(lon2d[lon] - dist_lon, 3), round(lat2d[lat] - dist_lat, 3)),  # bottom-left corner
                    (round(lon2d[lon] + dist_lon, 3), round(lat2d[lat] - dist_lat, 3)),  # bottom-right corner
                    (round(lon2d[lon] + dist_lon, 3), round(lat2d[lat] + dist_lat, 3)),  # top-right corner
                    (round(lon2d[lon] - dist_lon, 3), round(lat2d[lat] + dist_lat, 3)),  # top-left corner
                ])

                # Add the polygon to the corresponding list
                if mask[lat, lon] == 1:
                    urban_polygons.append(square)
                elif mask[lat, lon] == 0:
                    non_urban_polygons.append(square)
    else:
        dist_lat = abs(lat2d[1, 0] - lat2d[0, 0])/2
        dist_lon = abs(lon2d[0, 1] - lon2d[0, 0])/2

        # Iterate through the mask and generate polygons for urban (1) and non-urban (0) cells
        for lat in range(mask.shape[0] - 1):  # Avoid the last index to prevent out-of-bounds errors
            for lon in range(mask.shape[1] - 1):
                # Create a polygon using the 2D lat/lon coordinates of the cell corners
                if pd.isnull(mask.lon[lat, lon]): # If cell contains nans continue
                    continue
                square = Polygon([
                    (mask.lon[lat, lon] - dist_lon, mask.lat[lat, lon] - dist_lat),          # bottom-left corner
                    (mask.lon[lat, lon + 1] - dist_lon, mask.lat[lat, lon + 1] - dist_lat),  # bottom-right corner
                    (mask.lon[lat + 1, lon + 1] - dist_lon, mask.lat[lat + 1, lon + 1] - dist_lat),  # top-right corner
                    (mask.lon[lat + 1, lon] - dist_lon, mask.lat[lat + 1, lon] - dist_lat),  # top-left corner
                ])
                # Add the polygon to the corresponding list
                if mask[lat, lon] == 1:
                    urban_polygons.append(square)
                elif mask[lat, lon] == 0:
                    non_urban_polygons.append(square)
                    
    # Unite all adjacent polygons for urban (mask == 1) and non-urban (mask == 0)
    unified_urban_polygon = unary_union(urban_polygons)
    unified_non_urban_polygon = unary_union(non_urban_polygons)
    # Create GeoDataFrames for the urban and non-urban polygons
    # CRS 'EPSG:4326' specifies the WGS 84 coordinate system, which is widely used for global GPS coordinates (lat/lon)
    gdf_urban = gpd.GeoDataFrame(geometry=[unified_urban_polygon], crs='EPSG:4326')
    gdf_non_urban = gpd.GeoDataFrame(geometry=[unified_non_urban_polygon], crs='EPSG:4326')
    # Plot the boundary of the unified non-urban polygon
    # aspect='auto' prevents geopandas from overriding the axes aspect ratio
    # (which can fail when the geometry is empty or the axes is shared)
    if not gdf_non_urban.geometry.is_empty.all():
        gdf_non_urban.boundary.plot(ax=ax, color="#6B8D28", zorder=1, linewidth=2, aspect=None)
    # Plot the boundary of the unified urban polygon on top of the non-urban
    if not gdf_urban.geometry.is_empty.all():
        gdf_urban.boundary.plot(ax=ax, color="#9C3333", zorder=100, linewidth=2, aspect=None)

    return(gdf_urban, gdf_non_urban)


def plot_urban_borders(ds, ax, alpha = 1, linewidth = 2):
    """
    Plot the borders of urban areas on a map.

    Parameters:
    ----------
    ds (xr.Dataset): The dataset containing longitude, latitude, and urban area data.
    ax (matplotlib.axes._subplots.AxesSubplot): The matplotlib axes on which to plot.
    alpha (float): The transparency level of the borders (default is 1).
    linewidth (float): The width of the border lines (default is 2).

    Returns:
    -------
    None

    """
    lon2d = ds.lon.values
    lat2d = ds.lat.values
    dist_lat = lat2d[1, 0] - lat2d[0, 0]
    dist_lon = lon2d[0, 1] - lon2d[0, 0]
    dist_latlon = lat2d[0 ,1] - lat2d[0, 0]
    dist_lonlat = lon2d[1, 0] - lon2d[0, 0]
    # Overlay the cell borders and handle NaNs
    for i in range(len(ds.lat)-1):
        for j in range(len(ds.lon)-1):
            lons = [lon2d[i, j], lon2d[i, j+1], lon2d[i+1, j+1], lon2d[i+1, j], lon2d[i, j]]
            lats = [lat2d[i, j], lat2d[i, j+1], lat2d[i+1, j+1], lat2d[i+1, j], lat2d[i, j]]

            lons = lons - abs(lon2d[i, j] - lon2d[i, j+1])/2              
            lats = lats - abs(lat2d[i, j] - lat2d[i+1, j])/2
            
            data_cell = ds.values[i, j]
            
            if data_cell == 1:
                ax.plot(lons, lats, color='grey', zorder = 100, linewidth = linewidth, alpha = alpha)
            elif data_cell == 0:
                ax.plot(lons, lats, color='green', zorder = 1, linewidth = linewidth, alpha = alpha)

     # Plot the rightmost column
    for i in range(len(ds.lat) - 1):
        lons = [lon2d[i, -1], lon2d[i + 1, -1], lon2d[i + 1, -1] + dist_lon, lon2d[i, -1] + dist_lon, lon2d[i, -1]]
        lats = [lat2d[i, -1], lat2d[i + 1, -1], lat2d[i + 1, -1] + dist_latlon, lat2d[i, -1] + dist_latlon, lat2d[i, -1]]

        lons = lons - abs(lon2d[i, -1] - lon2d[i, -1])/2 - dist_lon/2            
        lats = lats - abs(lat2d[i, -1] - lat2d[i+1, -1])/2 
        
        data_cell = ds.values[i, -1]
    
        if data_cell == 1:
            ax.plot(lons, lats, color='grey', zorder=100, linewidth = linewidth, alpha = alpha)
        elif data_cell == 0:
            ax.plot(lons, lats, color='green', zorder=1, linewidth = linewidth, alpha = alpha)
    
    # Plot the topmost row
    for j in range(len(ds.lon) - 1):
        lons = [lon2d[-1, j], lon2d[-1, j + 1], lon2d[-1, j + 1]  + dist_lonlat, lon2d[-1, j]  + dist_lonlat, lon2d[-1, j]]
        lats = [lat2d[-1, j], lat2d[-1, j + 1], lat2d[-1, j + 1] + dist_lat, lat2d[-1, j] + dist_lat, lat2d[-1, j]]

        
        lons = lons - abs(lon2d[-1, j] - lon2d[-1, j+1])/2           
        lats = lats - abs(lat2d[-1, j] - lat2d[-1, j])/2 - dist_lat/2  
        
        data_cell = ds.values[-1, j]
    
        if data_cell == 1:
            ax.plot(lons, lats, color='red', zorder=100, linewidth=2)
        elif data_cell == 0:
            ax.plot(lons, lats, color='b', zorder=1, linewidth=2)
    
    # Plot the bottom right corner
    lons = [
        lon2d[-1, -1],
        lon2d[-1, -1] + dist_lon,
        lon2d[-1, -1] + dist_lon  + dist_lonlat,
        lon2d[-1, -1]  + dist_lonlat,
        lon2d[-1, -1]
    ]
    lats = [
        lat2d[-1, -1],
        lat2d[-1, -1] + dist_latlon,
        lat2d[-1, -1] + dist_lat + dist_latlon,
        lat2d[-1, -1] + dist_lat,
        lat2d[-1, -1]
    ]
        
    data_cell = ds.values[-1, -1]

    lons = lons - abs(lon2d[-1, -1] - lon2d[-1, -1])/2 - dist_lon/2
    lats = lats - abs(lat2d[-1, -1] - lat2d[-1, -1])/2 - dist_lat/2  
    
    if data_cell == 1:
        ax.plot(lons, lats, color="#6B8D28", zorder=100, linewidth=2)
    elif data_cell == 0:
        ax.plot(lons, lats, color="#9C3333", zorder=1, linewidth=2)

