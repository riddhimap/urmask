"""
Created by:
Riddhima Puri

This module contains the UrbanRuralSelection class, providing methods for selecting urban and, optionally, rural areas for high-resolution climate data.

It builds on the methodology developed by Diez-Sierra et al. (2025) in the URCLIMASK package (distributed under the Apache License 2.0), with modified logic 
for selecting urban and rural areas based on urban fraction, land fraction, and orography data.

It requires the following input datasets:
- Boundary dataset: GeoDataFrame containing city polygons (EPSG:4326) for the target city and neighboring cities.
- Urban fraction (uf) dataset: Contains the urban fraction values for each grid cell.
- Land fraction (lf) dataset: Contains the land fraction values for each grid cell.
- Orography (orog) dataset: Contains the orography values for each grid cell (optional, for rural mask creation).

References:
Diez-Sierra, J., Quintana, Y., Langendijk, G.S. et al. (2025).
A global CORDEX-based dataset delineating urban areas and their surroundings to assess climate change in megacities.
Scientific Data, 12, 1961.
https://doi.org/10.1038/s41597-025-06257-1

Diez-Sierra, J., Quintana, Y., Langendijk, G.S. et al. (2025). 
URCLIMASK: A Python Package for Delineating Urban Areas and Their Surrounding Reference Rural Regions from Regional Climate Models (RCMs) (Version v1.1.0) [Computer software]. 
Zenodo. https://doi.org/10.5281/zenodo.17257445
"""

import json

import matplotlib
import numpy as np
import xarray as xr
import geopandas as gpd

from matplotlib.colors import Colormap

import shapely.geometry
from shapely.prepared import prep
from shapely.geometry import Point,box
from skimage import measure
import skimage.morphology as morphology

from rapidfuzz import process


class UrbanRuralSelection:
    def __init__(
        self,
        *,
        city_name: str,
        all_cities_gdf: gpd.GeoDataFrame,
        gdf_column_city_names: str,
        all_cities_gdf_short_name: str | None = None,
        all_cities_gdf_long_name: str | None = None,

        resolution_km: float,
        model: str | None = None,
        domain: str | None = None,

        urban_threshold: float | None = 0.1,
        rural_threshold: float | None = 0.1,
        landsea_threshold: float = 70,
        orography_diff: float | None = 100,

    ):
        """
        Hyperparameters required for urban area selection.
    
        Parameters
        ----------
        city_name : str
            Name of the city to create mask for. Must match a city name in the GeoDataFrame.
        all_cities_gdf : gpd.GeoDataFrame
            GeoDataFrame containing all city polygons (EPSG:4326).
        gdf_column_city_names : str
            Name of the column containing city names in the GeoDataFrame.
        all_cities_gdf_short_name : str | None
            Short name or file name for the city polygons GeoDataFrame, used for metadata. Optional metadata.
        all_cities_gdf_long_name : str | None
            Long name or file name for the city polygons GeoDataFrame, used for metadata. Optional metadata.
        resolution_km : float
            Horizontal resolution of the dataset in kilometers. Used for buffering and distance calculations.
        model : str | None
            Name of the model or dataset being used. Optional metadata.
        domain : str | None
            Name of the domain or region being analyzed. Optional metadata.
        urban_threshold : float | None      
            Threshold of urban fraction in a grid cell above which it is considered urban. Default is >0.1.
            Must be same as low of threshold levels under create_urban_mask.
        rural_threshold : float | None
            Threshold of urban fraction in a grid cell below which it is considered rural. Default is <=0.1.
        landsea_threshold : float
            Threshold for percentage of land in a grid cell above which it is considered land (to exclude sea cells). Default is > 70%.
        orography_diff : float | None
            Orography difference threshold in meters to define rural areas. 
            If the elevation difference within the urban mask is less than this threshold, it will be used to define the rural mask. 
            Default is 100 m.
        """

        # Store parameters as instance variables

        # Check if target_city_name and city_name_column are provided and not empty
        if city_name is None or city_name.strip() == "":
            raise ValueError("city_name cannot be None or empty. Please provide a valid city name.")
        else:
            self.city_name = city_name
        
        if gdf_column_city_names is None or gdf_column_city_names.strip() == "":
            raise ValueError("gdf_column_city_names cannot be None or empty. Please provide a valid city name column.")
        else:
            self.gdf_column_city_names = gdf_column_city_names

        # Load city polygon from GeoDataFrame
        if isinstance(all_cities_gdf, gpd.GeoDataFrame):

            all_cities_gdf = all_cities_gdf.to_crs("EPSG:4326")  # ensure it's in geographic coordinates
            self.all_cities_gdf = all_cities_gdf # store the GeoDataFrame for future use

            # Get all unique city names from the specified column in the GeoDataFrame
            all_cities_names = (
                all_cities_gdf[gdf_column_city_names]
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )
            if not all_cities_names:
                raise ValueError(f"No valid city names found in column '{gdf_column_city_names}'.")

            # Do fuzzy matching to find city name in the GeoDataFrame
            match, score, _ = process.extractOne(city_name, all_cities_names)
            if score >= 80:  # threshold for fuzzy matching
                print(f"{city_name} exists in column '{gdf_column_city_names}' as '{match}'.")
                city_gdf = all_cities_gdf[(all_cities_gdf[gdf_column_city_names] == match)]
                self.city_name = match  # update city_name to the matched name
                if city_gdf.empty:
                    raise LookupError(f"No city polygon found for {match} in column '{gdf_column_city_names}'.")
                else:
                    self.city_polygon = city_gdf.geometry.union_all()  # merge all geometries if multiple entries exist
            else:
                raise LookupError(f"City '{city_name}' not found in column '{gdf_column_city_names}'.")
        else:
            raise TypeError("all_cities_gdf must be a GeoDataFrame.")
        

        self.all_cities_gdf_short_name = all_cities_gdf_short_name
        self.all_cities_gdf_long_name = all_cities_gdf_long_name

        self.resolution_km = resolution_km
        self.model = model
        self.domain = domain

        self.urban_th = urban_threshold
        self.rural_th = rural_threshold
        self.landsea_th = landsea_threshold
        self.orog_diff = orography_diff

        # initialize key variables for later use
        self.neighbors_gdf = None
        self.crop_distance_from_city_polygon_km = None
        self.rural_to_urban_ratio = None
        self.rural_area_inside_city_polygon = None
        self.urban_elev_min_threshold = None
        self.urban_elev_max_threshold = None
        self.thresholds = None


    def crop_city_area(self, *, 
                       ds: xr.DataArray,
                       crop_distance_from_city_polygon_km: float | None = None,
                       ) -> xr.DataArray:
        """
        Crop the area around a given city using city's polygon.

        Parameters:
        -----------
        ds : xarray.DataArray
            Dataset containing latitude and longitude coordinates.
        crop_distance_from_city_polygon_km : float
            Distance in kilometers to crop around the bounding box of the city polygon. 
            Default is approximately half the size of the bounding box.
        
        Returns:
        --------
        ds_cropped: xarray.DataArray
            Cropped dataset.
        """

        # Step 1: Check if city_polygon is set and if the dataset intersects with the city polygon
        if self.city_polygon is None:
            raise ValueError(
                "city_polygon is not set. Provide city_name or city_polygon when initializing UrbanSelection."
            )
        
        # check whether the dataset intersects with the city polygon 
        ds_bounds = box(
            float(ds.lon.min()), float(ds.lat.min()), 
            float(ds.lon.max()), float(ds.lat.max())
        )
        if not self.city_polygon.intersects(ds_bounds):
            raise ValueError(f"The dataset does not intersect with the city polygon for {self.city_name}. "
                             f"Please check the coordinates and the city polygon.")
                
        # Step 2: Find the bounding box around the city polygon
        lon_min, lat_min, lon_max, lat_max = self.city_polygon.bounds # (minx, miny, maxx, maxy)

        # Step 3: Crop the dataset to the bounding box of the city polygon plus an optional buffer distance in kilometers given by crop_distance_from_city_polygon_km. 
        if crop_distance_from_city_polygon_km is not None:
            self.crop_distance_from_city_polygon_km = crop_distance_from_city_polygon_km
            # Calculate buffer in terms of grid cells based on the resolution of the dataset
            buffer_cells = int(np.ceil(self.crop_distance_from_city_polygon_km / self.resolution_km)) # round or ceil
            print(f"Cropping dataset to city polygon with buffer of {self.crop_distance_from_city_polygon_km} km ({buffer_cells} grid cells) around the city polygon.")
        else:
            # If not provided, the buffer is set to half the size of the bounding box under Step 4.
            buffer_cells = None

        # Step 4: Crop dataset to the bounding box of the city polygon plus the buffer distance in grid cells. 
        # Handles both 1D and 2D lat/lon grids. Crops by indices rather than values to avoid issues with non-uniform grids or floating point precision.
        if ds["lat"].ndim == 1 and ds["lon"].ndim == 1:

            lat_vals = ds["lat"].values
            lon_vals = ds["lon"].values

            # Find indices inside polygon bounding box
            lat_idx = np.where((lat_vals >= lat_min) & (lat_vals <= lat_max))[0]
            lon_idx = np.where((lon_vals >= lon_min) & (lon_vals <= lon_max))[0]

            if len(lat_idx) == 0 or len(lon_idx) == 0:
                raise ValueError("No grid cells found inside bounds.")

            if buffer_cells is None:
                # Calculate number of grid cells covered by lat_idx and lon_idx
                num_lat_cells = lat_idx.max() - lat_idx.min() + 1
                num_lon_cells = lon_idx.max() - lon_idx.min() + 1
                buffer_cells = int(round(max(num_lat_cells, num_lon_cells) / 2))  # buffer half the size of the bounding box
                print(f"No 'crop_distance_from_city_polygon_km' provided. Cropping dataset to city polygon with buffer of {buffer_cells} grid cells around the city polygon.")
                self.crop_distance_from_city_polygon_km = buffer_cells * self.resolution_km

            # Expand indices (same logic as 2D case)
            i_min = max(0, lat_idx.min() - buffer_cells)
            i_max = min(len(lat_vals), lat_idx.max() + buffer_cells + 1)

            j_min = max(0, lon_idx.min() - buffer_cells)
            j_max = min(len(lon_vals), lon_idx.max() + buffer_cells + 1)

            # Slice the dataset 
            ds_cropped = ds.isel(lat=slice(i_min, i_max), lon=slice(j_min, j_max))

        elif ds["lat"].ndim == 2 and ds["lon"].ndim == 2:
            # Create a boolean mask for grid cells within the bounding box
            mask = (
                (ds["lon"] >= lon_min) & (ds["lon"] <= lon_max) &
                (ds["lat"] >= lat_min) & (ds["lat"] <= lat_max)
            )
            if not bool(mask.any()):
                raise ValueError("No grid cells found inside bounds.")
            
            # Find dimensions of lat and lon to slice the dataset
            dim_names = ds["lat"].dims

            # Get the row, column indices of the grid cells that fall within the bounding box
            indices = np.where(mask.values)

            if buffer_cells is None:
                # Calculate number of grid cells covered by the bounding box
                num_lat_cells = indices[0].max() - indices[0].min() + 1
                num_lon_cells = indices[1].max() - indices[1].min() + 1
                # print(f"Number of grid cells in bounding box: lat={num_lat_cells}, lon={num_lon_cells}")
                buffer_cells = int(round(max(num_lat_cells, num_lon_cells) / 2))  # buffer half the size of the bounding box
                print(f"No 'crop_distance_from_city_polygon_km' provided. Cropping dataset to city polygon with buffer of {buffer_cells} grid cells around the city polygon.")
                self.crop_distance_from_city_polygon_km = buffer_cells * self.resolution_km

            # Crop the dataset using the indices of the grid cells that fall within the bounding box
            # Compute bounds through indices (works for both increasing and decreasing)
            i_min = max(0, indices[0].min() - buffer_cells)
            i_max = min(ds.sizes[dim_names[0]], indices[0].max() + buffer_cells +1) # +1 to include the last index

            j_min = max(0, indices[1].min() - buffer_cells)
            j_max = min(ds.sizes[dim_names[1]], indices[1].max() + buffer_cells +1)

            # Slice
            ds_cropped = ds.isel({
                dim_names[0]: slice(i_min, i_max),
                dim_names[1]: slice(j_min, j_max),
            })
        else:
            raise ValueError("Unsupported lat/lon dimensionality.")
        
        
        # Also crop the all cities gdf if ds_cropped was evaluated
        
        lat = ds_cropped.lat.values
        lon = ds_cropped.lon.values

        # Find the bounding box of the cropped dataset to limit the search for neighboring cities to only those that could influence the grid cells in the cropped dataset.
        lon_min, lon_max = float(np.min(lon)), float(np.max(lon))
        lat_min, lat_max = float(np.min(lat)), float(np.max(lat))
        grid_box = box(lon_min, lat_min, lon_max, lat_max)

        # Find neighboring cities based on the actual grid footprint
        neighbors_gdf = self.all_cities_gdf[self.all_cities_gdf.geometry.intersects(grid_box)]
        # Exclude the target city itself from the neighboring cities
        neighbors_gdf = neighbors_gdf[neighbors_gdf[self.gdf_column_city_names] != self.city_name]
        self.neighbors_gdf = neighbors_gdf if not neighbors_gdf.empty else None

        return ds_cropped
    

    def create_urban_mask(self, *, 
                               ds_uf: xr.DataArray,
                               ds_lf: xr.DataArray,
                               ds_orog: xr.DataArray | None = None, # required for creating rural mask
                               thresholds_fraction: list[float] = None, # low to high tiers
                               debug_select_labels: bool = False,
                               orog_method_1: bool = True,
                               orog_method_2: bool = False,
                               ) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray | None, list[np.ndarray], list[np.ndarray]]:
        
        """
        Get city urban fraction, land fraction, and orography (optional) mask.

        Parameters
        ----------
        ds_uf : xarray DataArray
            Dataset containing the urban fraction and latitude/longitude grids.
        ds_lf : xarray DataArray
            Dataset containing the land-sea fraction and latitude/longitude grids.
        ds_orog : xarray DataArray, optional
            Dataset containing the orography and latitude/longitude grids. If not provided, the orography mask will be None.
        thresholds_fraction : [low to high] list of floats
            Threshold levels to define low/medium/high urban fraction tiers.
            If not provided, thresholds are dynamically computed based on the local distribution of urban fraction values above the given urban threshold.
        debug_select_labels : bool
            If True, debug information will be printed during the label selection process.

        Returns
        -------
        urban_mask : 2D boolean array
            Binary mask of the extracted urban footprint from urban fraction data AND land-sea fraction data.
        landsea_mask : 2D boolean array   
            Boolean mask with land-sea threshold applied on land-sea fraction data.
        orog_mask : 2D boolean array or None
            Boolean mask with orography threshold applied on orography data. None if orography data is not provided.
        labels : list of 2D int arrays
            List of labeled connected components for each density tier (low to high).
        labels_sel : list of 2D boolean arrays
            List of selected masks for each density tier after polygon intersection and fallbacks (low to high).
        """

        # Extract latitude, longitude, and urban fraction values
        lat = ds_uf['lat'].values
        lon = ds_uf['lon'].values

        # Standardize UF to [0, 1] and LF to [0, 100] before processing
        ds_uf, ds_lf = self._standardize_fraction_inputs(ds_uf=ds_uf, ds_lf=ds_lf)

        # Step 1: Compute dynamic thresholds if not provided
        if thresholds_fraction is None:
            thresholds_fraction = self._get_dynamic_thresholds(ds_uf=ds_uf)
            print(f"Dynamic thresholds for {self.city_name}: {thresholds_fraction}")
        else:
            if np.any(np.diff(thresholds_fraction) < 0):
                raise ValueError("thresholds_fraction must be sorted ascending.")
            if thresholds_fraction[0] < 0 or thresholds_fraction[-1] > 1:
                raise ValueError("thresholds_fraction must be within [0, 1].")

        self.thresholds = thresholds_fraction

        # Step 2: Compute masks for each urban fraction density level
        threshold_masks = self._compute_threshold_masks(
            ds_uf=ds_uf,
            thresholds_fraction=self.thresholds)
        if not threshold_masks:
            raise ValueError("No threshold masks were generated. Check thresholds_fraction.")

        # Step 3: Label connected regions independently for each density level
        labels = [measure.label(mask, connectivity=1) for mask in threshold_masks] #If connectivity=None, a full connectivity of input.ndim is used.
        
        # Step 4: Select labels for each density level that intersect the polygon
        labels_sel = self._select_labels_by_threshold_level(labels=labels, 
                                                            target_polygon=self.city_polygon, 
                                                            neighbors_polygons = self.neighbors_gdf.geometry if self.neighbors_gdf is not None else None,
                                                            lat=lat, lon=lon,
                                                            debug = debug_select_labels)

        # Step 5: Combine selected labels across all density tiers to create the final urban mask
        combined = np.logical_or.reduce(labels_sel)
        urban_mask = combined

        # Optional: Remove tiny isolated fragments
        # urban_mask = morphology.remove_small_objects(urban_mask, min_size=2, connectivity=2)

        # Step 6: Apply land fraction threshold to exclude sea cells
        landsea_mask = ds_lf > self.landsea_th   
        # Combine urban mask with land-sea masks
        urban_mask = urban_mask * landsea_mask

        if int(urban_mask.sum().item()) == 0:
            raise ValueError("Urban mask contains zero cells after land-sea filtering.")

        orog_mask = None

        # Optional: Create orography mask if orography data is provided
        if ds_orog is not None:
            """ Modified from: 
                Diez-Sierra et al. (2025). URCLIMASK: A Python Package for Delineating Urban Areas and Their Surrounding Reference Rural Regions from Regional Climate Models (RCMs) (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.17257445
            """
            # Calculate the maximum and minimum elevation values within the urban mask
            urban_elev_max = ds_orog.where(urban_mask).max().item()
            urban_elev_min = ds_orog.where(urban_mask).min().item()
            disliv=urban_elev_max-urban_elev_min
            # If the elevation difference is less than the user-defined threshold, use the user-defined threshold instead
            if disliv < self.orog_diff: # This condition includes those cases in which the urban areas is defined by a single point
                disliv = self.orog_diff
            # Define the upper and lower thresholds for the orography mask based on the urban elevation range and the user-defined threshold
            upper_thresh = urban_elev_max + disliv 
            lower_thresh = urban_elev_min - disliv 
            # Create the orography mask based on the defined thresholds
            orog_mask1 = ds_orog < upper_thresh 
            orog_mask2 = ds_orog > lower_thresh
            orog_mask = orog_mask1 & orog_mask2
                
            self.urban_elev_min_threshold = lower_thresh 
            self.urban_elev_max_threshold = upper_thresh


        urban_mask = urban_mask.astype(int).where(urban_mask.astype(int) == 1, np.nan)
        # Return urban mask as xarray and add attributes
        urban_mask = self._masks_to_dataset(
                                masks={"urban_mask": urban_mask},
                                variable_long_names={
                                    "urban_mask": "Urban mask: 1 for urban, NaN for the rest of the domain."
                                },
                                )

        return urban_mask, landsea_mask, orog_mask, labels, labels_sel
        
    
    def create_rural_mask(self, *, 
                            ds_uf: xr.DataArray,
                            urban_mask: xr.DataArray | xr.Dataset,
                            landsea_mask: xr.DataArray,
                            orog_mask: xr.DataArray | None = None,
                            rural_to_urban_ratio: float = 2.0,
                            rural_area_inside_city_polygon: bool = True
                            ) -> xr.Dataset:
            
            """
            Calculating rural mask through iterative growth of the urban mask into eligible rural areas.
            Note: The rural mask will grow only until number of non-urban cells exceed or the number of urban cells.

            Parameters:
            ----------
            ds_uf : xarray DataArray
                Dataset containing the urban fraction and latitude/longitude grids.
            urban_mask : xarray DataArray or xarray Dataset
                Boolean mask of the extracted urban footprint from urban fraction data.
            landsea_mask : xarray DataArray
                Boolean mask with land area fraction threshold applied from land fraction data.
            orog_mask : xarray DataArray | None
                Boolean mask with orography threshold applied from orography data. If not provided, the orography mask will be None.
            rural_to_urban_ratio : float
                Ratio of rural to urban cells. The rural mask will grow until the number of non-urban cells exceeds this ratio times the number of urban cells.
                Default is 2.0, meaning the rural mask can grow until it has twice as many cells as the urban mask.
                Always capped at the total number of eligible rural cells in the dataset.
            rural_area_inside_city_polygon : bool
                If True, rural mask can grow inside the city polygon. If False, rural mask will be restricted to areas outside the city polygon.
                Important when the city polygon is large and the user wants to avoid rural mask growth inside the city area OR when urban threshold > rural threshold.

            Returns:
            --------
            urban_mask : xarray Dataset
                Dataset containing the urban and rural masks with attributes.
            """
 
            # Code inspired by:
            # Diez-Sierra et al. (2025). URCLIMASK: A Python Package for Delineating Urban Areas and Their Surrounding Reference Rural Regions from Regional Climate Models (RCMs) (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.17257445

            # Step 1: Initialize variables

            # Standardize UF to [0, 1] and LF to [0, 100] before processing
            ds_uf, _ = self._standardize_fraction_inputs(ds_uf=ds_uf)
            
            # Define urban seed for rural growth
            urban_mask = urban_mask['urban_mask'] if isinstance(urban_mask, xr.Dataset) else urban_mask
            if bool(urban_mask.isnull().any().item()):
                urban_mask = urban_mask.fillna(0)
                print("Warning: Urban mask contains NaN values. Filling NaN with 0 for rural mask growth. Ensure that urban mask values are 1. Or else give urban_mask as boolean mask.")
            urban_bool = urban_mask.astype(bool)
            urban_cells = int(urban_bool.sum().item())
            if urban_cells == 0:
                raise ValueError("Urban mask contains no urban cells. Cannot define rural mask without urban seed.")

            # Number of cells added in the surroundings or accepted beyond the seed (this is what we want to control with the ratio and the distance threshold)
            non_urban_cells = 0

            # 4-neighbor growth then optional 8-neighbor fallback
            kernel1 = np.array([[0, 1, 0],
                                [1, 1, 1],
                                [0, 1, 0]])
            kernel2 = np.array([[1, 1, 1],
                                [1, 1, 1],
                                [1, 1, 1]])

            # Step 2: Define eligibility for rural growth

            # They must be below the urban fraction threshold and satisfy land mask and if provided, orography mask. 
            # If orography mask is not provided, only land mask is used.
            if orog_mask is None:
                valid_land = landsea_mask.astype(bool)
            else:
                valid_land = orog_mask.astype(bool) & landsea_mask.astype(bool)

            eligible = (ds_uf <= self.rural_th) & valid_land

            # Step 3: Define seed for rural growth: urban core + transition bridge + polygon mask (if rural_area_inside_city_polygon is False)

            # Optional: Enforce that selected rural cells are not inside the city polygon (exact polygon mask).
            self.rural_area_inside_city_polygon = rural_area_inside_city_polygon
            if self.rural_area_inside_city_polygon == False: 
                if ds_uf['lat'].ndim == 1 and ds_uf['lon'].ndim == 1:
                    lon2d, lat2d = np.meshgrid(ds_uf['lon'].values, ds_uf['lat'].values)
                    grid_dims = ('lat', 'lon')
                    grid_coords = {'lat': ds_uf['lat'], 'lon': ds_uf['lon']}
                elif ds_uf['lat'].ndim == 2 and ds_uf['lon'].ndim == 2:
                    lat2d = ds_uf['lat'].values
                    lon2d = ds_uf['lon'].values
                    grid_dims = ds_uf['lat'].dims
                    grid_coords = {grid_dims[0]: ds_uf[grid_dims[0]], grid_dims[1]: ds_uf[grid_dims[1]]}
                else:
                    raise ValueError('Unsupported lat/lon dimensionality in create_rural_mask')
                lon_flat = lon2d.ravel()
                lat_flat = lat2d.ravel()
                inside = np.array([self.city_polygon.contains(Point(lo, la)) for lo, la in zip(lon_flat, lat_flat)])
                polygon_mask = xr.DataArray(inside.reshape(lat2d.shape), dims=grid_dims, coords=grid_coords)
                # Apply dilation to the polygon mask to create a buffer around the city polygon, ensuring that rural cells do not grow into the city area. The buffer size is determined by the kernel2 (8-neighbor) structure.
                polygon_mask = xr.apply_ufunc(morphology.dilation, polygon_mask.astype(int), kwargs={"footprint": kernel2}).astype(bool)
                # Using polygon mask as eligibility condition is breaking connectivity (since polygon interior might have eligible cells)
                # eligible = (ds_uf <= self.rural_th) & valid_land & (~polygon_mask)
                # instead use it as seed
                print("Rural selection constraint: cells intersecting city polygon are excluded.")
            else:
                polygon_mask = xr.zeros_like(ds_uf, dtype=bool)

            # Handle transition zone between urban and rural thresholds
            if self.urban_th > self.rural_th:
                # the band between the urban and rural thresholds is not eligible for rural growth, so we can warn the user that the rural mask will be limited to a smaller area as it can break connectivty on side and limit growth.
                # print(f"Warning: Urban threshold ({self.urban_th}) is greater than rural threshold ({self.rural_th}). Rural mask will be limited to areas below the rural threshold.")
                
                # OR we can add the transition zone - which is neither urban seed nor eligible rural
                transition_zone = (ds_uf > self.rural_th) & (ds_uf <= self.urban_th) & valid_land
                # Urban cells + immediate surrounding ring --> it limits transition_bridge to only cells close to urban core, avoiding long-distance connections.
                near_urban = xr.apply_ufunc(morphology.dilation, urban_bool.astype(int), kwargs={"footprint": kernel2}).astype(bool)
                transition_bridge = transition_zone & near_urban
            else:
                transition_bridge = xr.zeros_like(urban_bool, dtype=bool)

            # transition_bridge = xr.zeros_like(urban_bool, dtype=bool)

            # Seed for growth = urban core + local transition bridge + polygon mask (if rural_area_inside_city_polygon is False)
            seed_for_growth = urban_bool | transition_bridge | polygon_mask
            # Currently dilated data is initialized as the seed for growth
            dilated_data = seed_for_growth.astype(bool)


            # Step 4: Define convergence conditions to stop the growth of the rural mask
            self.rural_to_urban_ratio = rural_to_urban_ratio

            # (i): The growth will stop when the number of non-urban cells exceeds the number of urban            
            eligible_total = int((eligible).sum().item())
            target_rural_cells = int(np.ceil(rural_to_urban_ratio * urban_cells))
            target_rural_cells = min(target_rural_cells, eligible_total)  # Ensure we don't exceed available eligible cells

            # (ii): The rural area can obviously not grow beyond the crop area distance threshold from the city polygon, so we can use that as a stopping condition.
            rural_growth_distance_from_polygon_km = self.crop_distance_from_city_polygon_km
            # maximum number of iterations to prevent infinite loops in case of non-convergence
            # based on the distance threshold and the dataset resolution, so it is adapted to the dataset and prevents overgrowing too much beyond the city polygon
            # maximum expansion radius in grid-cell steps
            counter_max = int(np.ceil(rural_growth_distance_from_polygon_km / self.resolution_km))  # Convert distance from city in kilometers to grid cells based on the dataset resolution (assuming square grid cells)
            counter = 0

            # Add a stagnation counter to stop the growth if no new non-urban cells are added for a certain number of iterations
            stagnation_steps = 0
            stagnation_limit = 3 

            # Step 5: Iteratively grow the region until the number of non-urban cells exceeds the number of urban cells 
            # or the maximum number of iterations is reached
            while counter < counter_max:
                
                prev_non_urban = non_urban_cells

                # First try: kernel1 (cross, 4-neighbor)
                candidate = xr.apply_ufunc(morphology.dilation, dilated_data, kwargs={"footprint": kernel1}).astype(bool)
                # Check eligibility and update candidate - keep urban seed always, and add eligible non-urban cells
                candidate = seed_for_growth | (candidate & eligible)
                # Count the number of non-urban cells in the candidate region
                new_non_urban = int((candidate & ~seed_for_growth).sum().item())

                # Fallback: if no growth and iterations remain, try kernel2 (square, 8-neighbor)
                if new_non_urban == prev_non_urban and counter < counter_max:
                    # Dilation with kernal2
                    candidate2 = xr.apply_ufunc(morphology.dilation, dilated_data, kwargs={"footprint": kernel2}).astype(bool)
                    # Check eligibility and update candidate2 - keep urban seed always, and add eligible non-urban cells
                    candidate2 = seed_for_growth | (candidate2 & eligible)
                    # Count number of non-urban cells in the candidate2 region
                    new_non_urban2 = int((candidate2 & ~seed_for_growth).sum().item())

                    if new_non_urban2 > new_non_urban:
                        candidate = candidate2
                        new_non_urban = new_non_urban2
                    elif new_non_urban2 == prev_non_urban:
                        print(f"Warning: No more non-urban cells found at iteration {counter}")
                        break

                dilated_data = candidate
                non_urban_cells = new_non_urban
                counter += 1
                print(f"Iteration {counter}: non-urban cells = {non_urban_cells}, urban cells = {urban_cells}")

                # success criterion
                if non_urban_cells >= target_rural_cells:
                    print("Stopping... Reached rural/urban ratio criterion.")
                    break

                # stagnation criterion
                if non_urban_cells == prev_non_urban:
                    stagnation_steps += 1
                else:
                    stagnation_steps = 0

                if stagnation_steps >= stagnation_limit:
                    print(f"Stopping... no growth for {stagnation_limit} iterations.")
                    break
            else:
                print("Stopping... reached iteration limit for rural growth.")
            
            # Step 6: Build final urban/vicinity mask
            # Non urban mask is the dilated data minus the seed for growth (urban + transition + polygon mask)
            non_urban_mask = (dilated_data & ~seed_for_growth).astype(int)
            # Assign 1 to urban, 0 to rural, and NaN to rest
            urban_rural_mask = urban_mask.astype(int).where(urban_mask.astype(int) == 1, np.nan)
            urban_rural_mask = urban_rural_mask.where(non_urban_mask.astype(int) == 0, 0)

            # Step 7: Return as xarray dataset with attributes
            return self._masks_to_dataset(
                masks={"urban_mask": urban_rural_mask},
                variable_long_names={
                    "urban_mask": "Urban mask: 1 urban, 0 rural vicinity, NaN rest of the domain."
                },
            )

    
    
    # ---- Helper functions -----

    def _standardize_fraction_inputs(
        self,
        *,
        ds_uf: xr.DataArray | None = None,
        ds_lf: xr.DataArray | None = None,
    ) -> tuple[xr.DataArray | None, xr.DataArray | None]:
        """
        Ensure urban fraction is in [0, 1] and land fraction is in [0, 100].
        Auto-converts when the data range indicates the wrong scale.

        Parameters
        ----------
        ds_uf : xr.DataArray | None
            Urban fraction data (expected 0–1 after standardization). Pass None to skip.
        ds_lf : xr.DataArray | None
            Land fraction data (expected 0–100 after standardization). Pass None to skip.

        Returns
        -------
        ds_uf, ds_lf : xr.DataArray | None
            Standardized copies of the inputs; None if not provided.
        """

        def _units_kind(da):
            u = (da.attrs.get("units", "") or "").strip().lower()
            if u in {"%", "percent", "percentage"}:
                return "percent"
            if u in {"fraction", "unitless", "dimensionless"}:
                return "fraction"
            return "unknown"

        if ds_uf is not None:
            uf_min = float(ds_uf.min(skipna=True).item())
            uf_max = float(ds_uf.max(skipna=True).item())
            if uf_min < 0 or uf_max > 100:
                raise ValueError(f"Urban fraction out of physical bounds: min={uf_min}, max={uf_max}")

            u_kind = _units_kind(ds_uf)
            if (u_kind == "percent" and uf_max > 1.0) or (u_kind == "unknown" and uf_max > 1.0):
                ds_uf = ds_uf / 100.0
                print(f"Urban fraction units detected as percent. Converting to fraction by dividing by 100.")


        if ds_lf is not None:
            lf_min = float(ds_lf.min(skipna=True).item())
            lf_max = float(ds_lf.max(skipna=True).item())
            if lf_min < 0 or lf_max > 100:
                raise ValueError(f"Land-sea fraction out of physical bounds: min={lf_min}, max={lf_max}")

            l_kind = _units_kind(ds_lf)
            if l_kind == "fraction":
                ds_lf = ds_lf * 100.0
                print(f"Land-sea fraction units detected as fraction. Converting to percent by multiplying by 100.")
            elif l_kind == "unknown":
                lf_p95 = float(ds_lf.quantile(0.95, skipna=True).item())
                if lf_p95 <= 1.0:
                    ds_lf = ds_lf * 100.0
                    print(f"Land-sea fraction units detected as unknown but 95th percentile <= 1. Converting to percent by multiplying by 100.")

        return ds_uf, ds_lf

    def _get_dynamic_thresholds(self, *, ds_uf: xr.DataArray):
        """
        Compute dynamic thresholds based on the local distribution of urban fraction values above the given urban threshold.

        Parameters:
        -----------
        ds_uf : xarray.DataArray
            Dataset containing the urban fraction (CPRCM) and latitude/longitude grids.

        Returns:
        --------
        thresholds : list of floats
            Threshold values corresponding to the given percentiles.
        """
        uf = ds_uf.values
        # Exclude NaN values (e.g., sea or masked cells) before any percentile calculation
        uf_valid = uf[np.isfinite(uf)]

        if uf_valid.size == 0:
            raise ValueError("No valid (non-NaN) urban fraction values found in ds_uf.")

        # Find where the urban threshold sits in the local valid-cell percentile distribution of the urban fraction values
        base_percentile = np.mean(uf_valid <= self.urban_th) * 100

        # Guard against degenerate case where almost all cells are below the threshold
        remaining_tail = 100.0 - base_percentile
        if remaining_tail < 1.0:
            print(f"Warning: Urban threshold {self.urban_th} is above the {base_percentile:.1f}th percentile. "
                  "Very few cells exceed it; returning a single threshold.")
            return [self.urban_th]

        # # minimum thresholds including the base threshold
        # n_thresholds = 2
        # n_thresholds = max(n_thresholds, round((100 - base_percentile) / 10))  # around 10% increments in the remaining tail
        # threshold_percentiles = np.linspace(base_percentile, 100.0, num=n_thresholds, endpoint=False)
        # print(f"Base percentile for urban threshold {self.urban_th}: {base_percentile:.2f}%")
        # print(f"Computed threshold percentiles: {threshold_percentiles}")
        # # Compute the actual the actual threshold values corresponding to the computed percentiles corres
        # thresholds = np.percentile(uf, threshold_percentiles)
        # # Ensure the first threshold is the physical base threshold
        # thresholds[0] = self.urban_th

        # Compute number of tiers: at least 2, roughly one per 10% of the remaining tail
        n_thresholds = max(2, round(remaining_tail / 10))
        threshold_percentiles = np.linspace(base_percentile, 100.0, num=n_thresholds, endpoint=False)

        print(f"Base percentile for urban threshold {self.urban_th}: {base_percentile:.2f}%")
        print(f"Computed threshold percentiles: {threshold_percentiles}")

        # Compute the actual threshold values from valid cells only corresponding to the computed percentiles
        thresholds = np.nanpercentile(uf_valid, threshold_percentiles)

        # Pin the first threshold to the physical base value (overrides percentile rounding)
        thresholds[0] = self.urban_th

        # Remove duplicate thresholds that would produce empty tiers in _compute_threshold_masks
        thresholds = np.unique(thresholds)
        
        return thresholds
    

    def _compute_threshold_masks(self, *, ds_uf: xr.DataArray, thresholds_fraction: list[float] = None):
        """
        Compute masks for each density level in the given thresholds for the urban fraction data.

        Parameters:
        -----------
        ds_uf : xarray.DataArray
            Dataset containing the urban fraction (CPRCM) and latitude/longitude grids.
        thresholds_fraction : list of floats, optional
            Threshold levels to define low/medium/high UF tiers. If not provided, default thresholds will be used.

        Returns:    
        --------    
        threshold_masks : list of xarray.DataArray 
            List of boolean masks for each threshold level.
        """
        uf = ds_uf.values
        threshold_masks = [(uf > thresholds_fraction[i]) & (uf <= thresholds_fraction[i + 1]) 
                           for i in range(len(thresholds_fraction) - 1)]
        if thresholds_fraction[-1] < 1:
            threshold_masks.append(uf > thresholds_fraction[-1])
        # print(threshold_masks)
        return threshold_masks
    
    
    def _select_labels_by_threshold_level(self, *, labels, target_polygon, neighbors_polygons, lat, lon, debug=False, debug_every=1):
        """
        For each threshold level, decide which connected-component labels belong to
        the target city polygon (and not to any neighbouring city polygon).

        Each label is first checked for a coarse intersection with the target / neighbour
        polygons under multiple scenarios:
        Scenario A: Label intersects target polygon but not any neighbour polygon -> select label
        Scenario B: Label intersects neighbour polygon but not the target polygon -> discard label
        Scenario C: Label intersects both target and neighbour polygons -> resolve on a cell-by-cell basis using distance criteria:
                    1. Physical containment  (covers / within tolerance)
                    2. Boundary distance     (normalised to local grid steps)
                    3. Centroid distance     (tie-breaker)
                    4. Centroid x-coordinate (last-resort tie-breaker)

        Parameters
        ----------
        labels : list of 2-D int arrays
            Connected-component label arrays, one per threshold level.
        target_polygon : shapely.geometry.Polygon
            The city polygon whose urban cells we want to select.
        neighbors_polygons : list of shapely.geometry.Polygon or None
            Adjacent city polygons that may share connected components with the target.
        lat, lon : 1-D or 2-D arrays
            Grid coordinates.
        debug : bool
            Print per-label diagnostics when True.
        debug_every : int
            Print diagnostics only every Nth threshold level (reduces noise).

        Returns
        -------
        labels_sel : list of 2-D bool arrays
            Per-level boolean masks of selected cells.
        """

        # --- Initialize selection masks for each threshold level ---
        labels_sel = [np.zeros_like(lbl, dtype=bool) for lbl in labels]

        # --- Ensure 2D coordinates ---
        if lat.ndim == 1 and lon.ndim == 1:
            lon2d, lat2d = np.meshgrid(lon, lat)
        elif lat.ndim == 2 and lon.ndim == 2:
            lat2d, lon2d = lat, lon
        else:
            raise ValueError("Unsupported lat/lon array dimensions")

        # --- Local grid spacing (degrees) ---
        dlat_i, dlat_j = np.gradient(lat2d)
        dlon_i, dlon_j = np.gradient(lon2d)
        # Grid size = max gradient magnitude in both directions (degrees).
        grid_size = np.sqrt(
            np.maximum(np.abs(dlat_i), np.abs(dlat_j)) ** 2
            + np.maximum(np.abs(dlon_i), np.abs(dlon_j)) ** 2
        )
        # Use half the grid size as a tolerance for distance checks
        tol_grid = 1/2 * grid_size
        # Small epsilon to avoid numerical issues in distance comparisons
        eps = 1e-12

        # --- Precompute prepared geometries for efficient intersection tests ---
        prep_target = prep(target_polygon)
        target_center = target_polygon.centroid

        # Handle neighbors polygons if none are present
        if neighbors_polygons is None or len(neighbors_polygons) == 0:
            neighbors_polygons = []
            prep_neighbors = []
            neighbor_centers = []
        else:
            neighbors_polygons = list(neighbors_polygons)
            prep_neighbors = [prep(p) for p in neighbors_polygons]
            neighbor_centers = [p.centroid for p in neighbors_polygons]

        if debug:
            print(f"[DEBUG] n_levels={len(labels)}, n_neighbors={len(neighbors_polygons)}")


        # --- Iterate over each threshold level and evaluate labels against target and neighbor polygons ---
        for k, lbl in enumerate(labels):

            # Extract unique label IDs (skip the background label 0)
            ids = np.unique(lbl)
            ids = ids[ids != 0]

            if debug and (k % debug_every == 0):
                print(f"\n[DEBUG] level={k}, n_labels={len(ids)}")

            # Iterate over each unique label ID
            for lab in ids:

                # Create a mask for the current label
                mask = (lbl == lab)
                idx = np.where(mask)
                n_cells = idx[0].size

                # --- Intersection checks: whether any cell of the label intersects with the target polygon or any neighbor polygon ---
                intersects_target = False
                intersects_neighbor = False

                for i, j in zip(idx[0], idx[1]):
                    point = Point(lon2d[i, j], lat2d[i, j])

                    # Check if the point intersects with the target polygon or is within the local tolerance distance
                    if not intersects_target:
                        if prep_target.intersects(point) or (target_polygon.distance(point) <= tol_grid[i, j]):
                            intersects_target = True

                    # Check if the point intersects with any neighbor polygon or is within the local tolerance distance
                    if not intersects_neighbor:
                        for prep_n, pol_n in zip(prep_neighbors, neighbors_polygons):
                            if prep_n.intersects(point) or (pol_n.distance(point) <= tol_grid[i, j]):
                                intersects_neighbor = True
                                break

                    # If both intersections are found, no need to check further cells for this label
                    if intersects_target and intersects_neighbor:
                        break

                if debug and (k % debug_every == 0):
                    print(
                        f"[DEBUG] level={k}, lab={lab}, n_cells={n_cells}, "
                        f"intersects_target={intersects_target}, intersects_neighbor={intersects_neighbor}"
                    )

                # Label selection logic based on intersection scenarios:
                # Scenario A: Label intersects target polygon but not any neighbor polygon
                if intersects_target and not intersects_neighbor:
                    labels_sel[k][mask] = True
                    continue

                # Scenario B: Label intersects neighbor polygon but not the target polygon
                if not intersects_target and intersects_neighbor:
                    continue

                # No neighbors at all - if the label intersects the target polygon, it is directly selected
                if intersects_target and (len(neighbors_polygons) == 0):
                    labels_sel[k][mask] = True
                    continue

                # Scenario C: Label intersects BOTH target and neighbor polygons -> requires further checking
                if intersects_target and intersects_neighbor:
                    assigned_true = 0
                    assigned_false = 0
                    ambiguous = 0

                    # Check whether each cell of the label is inside the target polygon or any neighbor polygon, and evaluate distances to determine assignment
                    for i, j in zip(idx[0], idx[1]):
                        point = Point(lon2d[i, j], lat2d[i, j])

                        # Check if the point is inside the target polygon and its distance to the target polygon
                        in_target = prep_target.covers(point)
                        dist_to_target = target_polygon.distance(point)

                        # Check if the point is inside any neighbor polygon and find the closest neighbor polygon and its distance
                        min_neighbor_dist = np.inf
                        closest_neighbor_idx = -1
                        in_neighbor = False

                        for n_idx, pol_n in enumerate(neighbors_polygons):
                            d = pol_n.distance(point)
                            if d < min_neighbor_dist:
                                min_neighbor_dist = d
                                closest_neighbor_idx = n_idx
                            if prep_neighbors[n_idx].covers(point):
                                in_neighbor = True


                        # Step 1: cell is inside (or on boundary of) the target
                        # A cell physically inside the target is never given away.
                        # The only exception: cell is inside the neighbour AND outside the target AND the neighbour boundary is closer -> yield to neighbour.
                        if in_target or (dist_to_target <= tol_grid[i, j]): # similar as the check for label intersection!

                            # Only yield to neighbor if neighbor is a significantly better fit for it
                            # not in_target ensures that if a grid cell is physically inside the target city, it cannot be given away to a neighbor based on distance.
                            if in_neighbor and not in_target and (min_neighbor_dist < dist_to_target):
                                if debug:
                                    print(f"[DEBUG] level={k}, lab={lab}, n_cells={n_cells} - covers_target={in_target}, covers_neighbor={in_neighbor} -> assigned False (neighbor closer)")
                                labels_sel[k][i, j] = False
                                assigned_false += 1
                        
                            else:
                                if debug:
                                    print(f"[DEBUG] level={k}, lab={lab}, n_cells={n_cells} - covers_target={in_target}, covers_neighbor={in_neighbor} -> assigned True (target closer or equal)")
                                labels_sel[k][i, j] = True
                                assigned_true += 1
                            continue

                        # Step 2: cell it outside both polygons
                        # Normalise distances by the local grid step so the comparison is
                        # scale-independent.  Assign to whichever polygon boundary is closer.

                        # Distance in local grid steps
                        local_cell_size = max(grid_size[i, j], eps)
                        steps_to_target = target_polygon.distance(point) / local_cell_size
                        steps_to_neighbor = min_neighbor_dist / local_cell_size

                        # If cell is closer to target than to the neighbor, assign to taget
                        if steps_to_target < steps_to_neighbor:
                            labels_sel[k][i, j] = True
                            assigned_true += 1

                        # Step 3: If cell is equidistance to both target and neighbor, need to break the tie 
                        # Can use the distance to the centroids of the polygons as a criterion.
                        elif np.isclose(steps_to_target, steps_to_neighbor, rtol=0.0, atol=1e-9):
                            ambiguous += 1
                            # If closest_neighbor_idx is -1, it means there are no neighbors, so we should assign to target by default
                            if closest_neighbor_idx < 0:
                                in_target_bool = in_target or (dist_to_target <= tol_grid[i, j])
                                labels_sel[k][i, j] = in_target_bool
                            else:

                                cos_lat = np.cos(np.radians(lat2d[i, j]))
                                dlon_target = (point.x - target_center.x) * cos_lat
                                dlat_target = (point.y - target_center.y)
                                dist_target = np.hypot(dlon_target, dlat_target)

                                n_center = neighbor_centers[closest_neighbor_idx]
                                dlon_neighbor = (point.x - n_center.x) * cos_lat
                                dlat_neighbor = (point.y - n_center.y)
                                dist_neighbor = np.hypot(dlon_neighbor, dlat_neighbor)

                                # if the distance to the target centroid is less than the distance to the neighbor centroid, assign to target
                                if dist_target < dist_neighbor:
                                    labels_sel[k][i, j] = True
                                # Step 4: if equidistant to both centroids, assign based on the x-coordinate (longitude) of the centroids as a tie-breaker
                                elif np.isclose(dist_target, dist_neighbor, rtol=0.0, atol=1e-9):
                                    labels_sel[k][i, j] = (target_center.x > n_center.x)
                                # if the distance to the neighbor centroid is less than the distance to the target centroid, assign to neighbor
                                else:
                                    labels_sel[k][i, j] = False
                                
                                # Debugging output for ambiguous cases
                                if labels_sel[k][i, j]:
                                    assigned_true += 1
                                else:
                                    assigned_false += 1
                        else:
                            # If the cell is closer to the neighbor than to the target, assign to neighbor
                            labels_sel[k][i, j] = False
                            assigned_false += 1

                    if debug and (k % debug_every == 0):
                        kept = int(labels_sel[k][mask].sum())
                        print(
                            f"[DEBUG] level={k}, lab={lab}, kept={kept}/{n_cells}, "
                            f"assigned_true={assigned_true}, assigned_false={assigned_false}, ambiguous={ambiguous}"
                        )

        return labels_sel
    

    def _masks_to_dataset(
        self,
        *,
        masks: dict[str, xr.DataArray | None],
        variable_long_names: dict[str, str] | None = None,
    ) -> xr.Dataset:
        """
        Convert a dictionary of masks into an xarray.Dataset and attach metadata.

        Parameters
        ----------
        masks : dict[str, xarray.DataArray | None]
            Mapping of variable name to mask DataArray. Entries with None are skipped.
        variable_long_names : dict[str, str], optional
            Optional mapping of variable name to long_name metadata.

        Returns
        -------
        xarray.Dataset
            Dataset containing provided masks and metadata attributes.
        """
        urmask = xr.Dataset()
        for var_name, mask in masks.items():
            if mask is None:
                continue
            # If already a DataArray, keep coordinates/dimensions unchanged.
            if isinstance(mask, xr.DataArray):
                urmask[var_name] = mask
            else:
                urmask[var_name] = mask.astype(int) if getattr(mask, "dtype", None) == bool else mask

        urmask = self._netcdf_attrs(urmask, variable_long_names=variable_long_names)
        return urmask


    def _netcdf_attrs(self, urmask, variable_long_names: dict[str, str] | None = None):
        """
        Add attributes to the output netCDF dataset for reproducibility and clarity.

        Parameters
        ----------
        urmask : xarray.Dataset
            The dataset containing the urban mask and related variables.
        variable_long_names : dict[str, str], optional
            Optional mapping of variable names to their long_name attributes.

        Returns
        -------
        urmask : xarray.Dataset
            The same dataset with added attributes.
        """

        urmask = urmask.copy()

        # Helper: convert values to netCDF-safe attribute types
        def _to_netcdf_attr_value(value):

            # Convert bool to string for widest netCDF compatibility
            if isinstance(value, bool):
                return str(value)
            
            # Keep plain Python scalars
            if isinstance(value, (str, int, float)):
                return value

            # numpy scalar -> Python scalar
            if isinstance(value, np.generic):
                py_val = value.item()
                return str(py_val) if isinstance(py_val, bool) else py_val

            # Lists/tuples/ndarrays -> JSON string
            if isinstance(value, (list, tuple, np.ndarray)):
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                return json.dumps(value)

            # Fallback for unknown objects
            return str(value)

        # Fixed dataset-level metadata
        urmask.attrs["created_by"] = "UrbanRuralSelection class in urmask.py"
        urmask.attrs["description"] = (
            "Dataset containing urban mask for the specified city. "
            "Mask values follow example convention (e.g., 1 urban, 0 rural vicinity, NaN outside domain)."
        )

        # Optional metadata from class attributes (only if present and not None) for reproducibility
        attrs_map = {
            "city_name": "city_name",
            "all_cities_gdf_short_name": "all_cities_gdf_short_name",
            "all_cities_gdf_long_name": "all_cities_gdf_long_name",
            "model": "model",
            "domain": "domain",
            "resolution_km": "resolution_km",
            "urban_area_urban_fraction_threshold": "urban_th",
            "urban_area_threshold_levels": "thresholds",
            "rural_area_urban_fraction_threshold": "rural_th",
            "landsea_fraction_threshold": "landsea_th",
            "orography_difference_threshold": "orog_diff",
            "crop_city_area_distance_from_city_polygon_km": "crop_distance_from_city_polygon_km",
            "rural_to_urban_ratio": "rural_to_urban_ratio",
            "rural_area_inside_city_polygon": "rural_area_inside_city_polygon",
        }
        
        for out_key, src_attr in attrs_map.items():
            value = getattr(self, src_attr, None)
            if value is not None:
                urmask.attrs[out_key] = _to_netcdf_attr_value(value)


        # Variable-level metadata
        default_long_names = {
            "urban_mask": "Urban-rural mask: 1 urban, 0 rural vicinity, NaN rest of the domain.",
            "landsea_mask": "Land-area mask: 1 land, 0 sea.",
            "orog_mask": "Orography mask: 1 indicates areas meeting orography criteria, 0 otherwise.",
        }
        if variable_long_names is not None:
            if not isinstance(variable_long_names, dict):
                raise ValueError("variable_long_names must be a dictionary mapping variable names to long_name strings.")
            default_long_names.update(variable_long_names)

            # Optional warning for unknown keys passed by caller
            unknown = set(variable_long_names.keys()) - set(urmask.data_vars.keys())
            if unknown:
                print(
                    f"Warning: variable_long_names contains keys not present in dataset: {sorted(unknown)}"
                )

        for var_name in urmask.data_vars:
            if var_name in default_long_names:
                urmask[var_name].attrs["long_name"] = default_long_names[var_name]

        # Convert boolean attributes to strings for netCDF4 compatibility
        for attr_name in list(urmask.attrs.keys()):
                urmask.attrs[attr_name] = _to_netcdf_attr_value(urmask.attrs[attr_name])

        return urmask


    # ---- Plotting functions -----

    def plot_variables_masks(self,
                            *,
                            ds_uf: xr.DataArray,
                            ds_lf: xr.DataArray,
                            ds_orog: xr.DataArray = None,
                            urban_mask: xr.DataArray | xr.Dataset = None,
                            landsea_mask: xr.DataArray = None,
                            orog_mask: xr.DataArray = None,
                            neighbors = False,
                            urban_cmap: str = "Grays",
                            landsea_cmap: str = "YlGnBu_r",
                            orog_cmap: str = None,
                            polygon_color: str = "#ff0000",
                            save_path: str = None,
                            ):
        """
        Plot the urban fraction, land-sea fraction, and orography datasets along with their corresponding masks and the city polygon boundary.

        Parameters
        ----------
        ds_uf : xarray.DataArray
            Dataset containing the urban fraction and latitude/longitude grids.
        ds_lf : xarray.DataArray
            Dataset containing the land-sea fraction and latitude/longitude grids.
        ds_orog : xarray.DataArray, optional
            Dataset containing the orography and latitude/longitude grids.
        urban_mask : xarray.DataArray or xarray.Dataset, optional
            Binary mask indicating urban/rural areas.
        landsea_mask : xarray.DataArray
            Binary mask indicating meeting land-sea areas criteria.
        orog_mask : xarray.DataArray, optional
            Binary mask indicating areas meeting orography criteria.
        urban_cmap : str, optional
            Colormap for urban fraction.
        landsea_cmap : str, optional
            Colormap for land-sea fraction.
        orog_cmap : str, optional
            Colormap for orography.
        polygon_color : str, optional
            Color for the city polygon boundary.
        save_path : str, optional
            Path to save the plot. If None, the plot will be displayed instead of saved.

        """
        # Lazy imports: plotting stack is heavy and should not slow module import.
        import matplotlib
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        from matplotlib.colors import LinearSegmentedColormap, ListedColormap, BoundaryNorm
        plt.rcParams.update({'axes.titlesize': 16, 'axes.labelsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14})

        # Standardize UF to [0, 1] and LF to [0, 100] before processing
        ds_uf, ds_lf = self._standardize_fraction_inputs(ds_uf=ds_uf, ds_lf=ds_lf)

        # Define colormaps and normalization for urban fraction, land-sea fraction, and orography
        base_cmap = self._resolve_cmap(urban_cmap, matplotlib_mod=matplotlib, name_for_error="urban_cmap")
        colors = base_cmap([0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        cmap_uf = ListedColormap(colors)
        boundaries = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        norm_uf = BoundaryNorm(boundaries, cmap_uf.N, extend='neither')

        base_cmap = self._resolve_cmap(landsea_cmap, matplotlib_mod=matplotlib, name_for_error="landsea_cmap")
        colors = base_cmap(np.linspace(0, 1, 11))
        cmap_lf = ListedColormap(colors)
        boundaries = np.arange(0, 110, 10)
        norm_lf = BoundaryNorm(boundaries, cmap_lf.N, extend='neither')

        if ds_orog is not None:
            if orog_cmap is None:
                colors = ['#278908', '#faf998', '#66473b']
                base_cmap = LinearSegmentedColormap.from_list("custom_terrain", colors)
            # elif isinstance(orog_cmap, Colormap):
            #     base_orog = orog_cmap
            # else:
            #     base_orog = matplotlib.colormaps.get_cmap(orog_cmap)
            else:
                base_cmap = self._resolve_cmap(orog_cmap, matplotlib_mod=matplotlib, name_for_error="orog_cmap")

            cmap_orog = ListedColormap(base_cmap(np.linspace(0, 1, 11)))
            orog_min = float(np.nanmin(ds_orog))
            orog_max = float(np.nanmax(ds_orog))
            orog_bounds = np.linspace(orog_min, orog_max, 11)
            norm_orog = BoundaryNorm(orog_bounds, cmap_orog.N, extend="neither")
    
        # Set up the figure and axes for plotting
        proj = ccrs.PlateCarree()

        if ds_orog is None and orog_mask is None:
            nrows, ncols = 2, 2
            figsize=(15, 12)
        else:
            nrows, ncols = 2, 3
            figsize=(20, 10)

        fig, axes = plt.subplots(nrows, ncols, subplot_kw={'projection': proj}, figsize=figsize)

        # Plot the city polygon on all subplots
        for l in range(nrows):
            for k in range(ncols):
                # self.city_polygon.boundary.plot(ax=axes[l, k], facecolor="none", transform=proj, edgecolor="#ff0000", linewidth=1, zorder=1000, label=f"{self.city_polygon_short_name}")
                axes[l, k].add_geometries([self.city_polygon], crs=proj, facecolor="none", edgecolor= polygon_color, linewidth=1, zorder=1000, label=f"{self.all_cities_gdf_short_name} {self.city_name}")
                if neighbors and self.neighbors_gdf is not None and len(self.neighbors_gdf) > 0:
                    neighbors_polygons = self.neighbors_gdf.geometry
                    axes[l, k].add_geometries(neighbors_polygons, crs=proj, facecolor="none", edgecolor=polygon_color, linestyle='--', linewidth=0.6, zorder=1000, label=f"{self.all_cities_gdf_short_name} Neighbors")
                axes[l, k].set_xlabel('Longitude')
                axes[l, k].set_ylabel('Latitude')
                axes[l, k].coastlines()
                # axes[l, k].gridlines(draw_labels=True)
    
        # Plot the urban fraction, land-sea fraction, and orography datasets on the first row of subplots
        im1 = axes[0, 0].pcolormesh(ds_uf.lon, ds_uf.lat,
                                    ds_uf,
                                    cmap=cmap_uf,
                                    norm=norm_uf,
                                    transform=proj)
        fig.colorbar(im1, ax=axes[0, 0], shrink=0.8, label='Urban Fraction')
        # axes[0, 0].set_title('Urban Fraction')


        im2 = axes[0, 1].pcolormesh(ds_lf.lon, ds_lf.lat,
                                    ds_lf,
                                    cmap=cmap_lf,
                                    norm=norm_lf,
                                    transform=proj)
        fig.colorbar(im2, ax=axes[0, 1], shrink=0.8, label='Land-Sea Fraction [%]')
        # axes[0, 1].set_title('Land-Sea Fraction [%]')


        if ds_orog is not None:
            im3 = axes[0, 2].pcolormesh(ds_orog.lon, ds_orog.lat,
                                    ds_orog,
                                    cmap=cmap_orog,
                                    norm=norm_orog,
                                    transform=proj)
            fig.colorbar(im3, ax=axes[0, 2], shrink=0.8, label='Orography [m]')
            # axes[0, 2].set_title('Orography [m]')

        # Plot the urban mask, land-sea mask, and orography mask on the second row of subplots
        if urban_mask is not None:
            urban_mask = urban_mask['urban_mask'] if isinstance(urban_mask, xr.Dataset) else urban_mask
            im4 = axes[1, 0].pcolormesh(ds_uf.lon, ds_uf.lat,
                                        ds_uf.where(urban_mask == 1, np.nan),
                                        cmap=cmap_uf,
                                        norm=norm_uf,
                                        transform=proj)
            fig.colorbar(im4, ax=axes[1, 0], shrink=0.8, label='Urban Fraction (UF)')
            axes[1, 0].set_title(f"Urban-Rural Mask \nUrban: > {self.urban_th} & Rural: <= {self.rural_th}") #, urb_sur_th = {self.urban_sur_th}\nratio_r2u = {self.ratio_r2u}, max_city = {self.min_city_size})")
            # axes[1, 0].set_title(f"Urban Mask: Urban Fraction > {self.urban_th} \n& Rural Mask: Urban Fraction <= {self.rural_th}") #, urb_sur_th = {self.urban_sur_th}\nratio_r2u = {self.ratio_r2u}, max_city = {self.min_city_size})")
        else:
            # Hide
            axes[1,0].set_visible(False)

        if landsea_mask is not None:
            im5 = axes[1, 1].pcolormesh(ds_lf.lon, ds_lf.lat,
                                        ds_lf.where(landsea_mask == 1, np.nan),
                                        cmap=cmap_lf,
                                        norm=norm_lf,
                                        transform=proj)
            fig.colorbar(im5, ax=axes[1, 1], shrink=0.8, label='Land-Sea Fraction [%]')
            axes[1, 1].set_title(f'Urban-Rural Mask \n Land-Sea Fraction > {self.landsea_th}%')
        else:
            # Hide
            axes[1,1].set_visible(False)

        if orog_mask is not None:
            im6 = axes[1, 2].pcolormesh(ds_orog.lon, ds_orog.lat,
                                        ds_orog.where(orog_mask == 1, np.nan),
                                        cmap=cmap_orog,
                                        norm=norm_orog,
                                        transform=proj)
            fig.colorbar(im6, ax=axes[1, 2], shrink=0.8, label='Orography [m]')
            elev_lim_min = self.urban_elev_min_threshold
            elev_lim_max = self.urban_elev_max_threshold
            axes[1, 2].set_title(f'Urban-Rural Mask \n Orography Difference = {self.orog_diff}m \n{elev_lim_min:.0f}m < Orography < {elev_lim_max:.0f}m')
        elif ds_orog is not None and orog_mask is None:
            # Hide
            axes[1,2].set_visible(False)

        # plot the urban and rural masks as polygons on top for better visualization
        if urban_mask is not None:
            urban_mask = urban_mask['urban_mask'] if isinstance(urban_mask, xr.Dataset) else urban_mask
            from .utils import plot_urban_polygon
            for k in range(ncols):
                plot_urban_polygon(urban_mask, axes[1, k])
        

        # plt.subplots_adjust(wspace=0.1, hspace=0.1)
        # Only add legend if there are labeled artists (suppress warning for CartoPy features)
        handles, labels = plt.gca().get_legend_handles_labels()
        if labels:
            plt.legend(
                loc='lower center',  # Place the legend at the lower center
                bbox_to_anchor=(0.5, -0.2),  # Adjust the position to be outside the figure
                ncol=3,  # Arrange the legend items in a single column
                fontsize=10,
            )
        fig.suptitle(f'{self.city_name}', fontsize=18, y=0.9, fontweight='bold')  # Adjust y for better placement

        # Save the figure if a save path is provided
        if save_path is not None:
            plt.savefig(save_path, bbox_inches='tight', dpi=300)
            print(f"Figure saved to {save_path}")            
        return fig
        
# ----- Helper function to resolve colormap input -----
    def _resolve_cmap(self, cmap_input, *, matplotlib_mod, name_for_error: str):
        if isinstance(cmap_input, Colormap):
            return cmap_input
        if isinstance(cmap_input, str):
            if cmap_input in matplotlib_mod.colormaps:
                return matplotlib_mod.colormaps.get_cmap(cmap_input)
            raise ValueError(f"Unknown colormap name: {cmap_input}")
        raise TypeError(f"{name_for_error} must be a colormap name or a Colormap object")
