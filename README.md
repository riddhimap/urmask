# urmask: Urban and Rural Area Classification for High-Resolution Climate Model Data

This repository contains the _**urmask**_ Python package, designed to derive urban and, optionally, rural masks from climate model data. 

_urmask_ builds upon the methodology introduced in URCLIMASK (Diez-Sierra et al., 2025) and extends it for convection-permitting regional climate models (CPRCMs) operating at kilometre-scale resolutions.

URCLIMASK was developed primarily for coarse-resolution regional climate models (12.5–25 km), where cities are represented by only a few grid cells. At CPRCM resolutions (~3 km), urban areas become much more spatially complex, with urban cores, suburbs, satellite towns, and neighboring cities often being explicitly resolved. Thus, at this scale, nearby urban areas can appear connected, making it difficult to identify independent cities with a simple single-threshold approach.

To address these challenges, _urmask_ introduces:

- Multi-threshold urban-density analysis to improve the separation of connected urban areas.
- Boundary-informed urban area extraction to distinguish the target city from neighboring cities.
- Automated parameter selection and convergence criteria that reduce manual parameter tuning.
- Unique grid-cell assignment, ensuring that a grid cell cannot belong to multiple cities and enabling robust inter-city and regional comparisons.

Unlike a simple administrative clipping approach, _urmask_ retains the land-based approach of URCLIMASK by using the model's urban fraction, land-sea fraction, and orography while incorporating administrative boundaries as additional spatial context. In doing so, _urmask_ extends the original methodology for the greater spatial complexity of kilometre-scale climate simulations.

The framework supports climate-service applications and can be adapted to different use cases. For example, users may define urban areas using City (C) or Functional Urban Areas (F) boundaries from the Urban Audit dataset (URAU) (European Commission, Eurostat, 2025). This approach enables the generation of administratively meaningful urban and rural masks that remain consistent with the land-cover representation of individual climate models.

A key advantage of _urmask_ is that many parameters can be derived dynamically from the characteristics of the target city and the underlying model data. This minimizes the need for city-specific parameter tuning, allowing the methodology to be applied consistently across cities of different sizes and morphologies while improving robustness, transferability, and comparability.

The result is thus a transferable workflow for generating urban and rural masks that is particularly suited to high-resolution climate simulations, Urban Heat Island studies, climate service applications, and inter-city climate analyses.

## Repository Contents

### **urban_rural_mask_example.ipynb**

A Jupyter Notebook demonstrating the complete workflow for generating urban and rural masks using:

- A regional climate model land-representation dataset.
- A city or urban boundary polygon.

_Urban fraction, land-sea fraction, and orography datasets are not distributed with this repository. Examples shown are for the RCMs used as part of the NUKLEUS simulations (Sieck et al., 2026). Dataset available upon request._

_Urban boundary datasets are not distributed with this repository and must be downloaded separately from the sources listed below._


### **urmask/urmask.py**

Contains the UrbanRuralSelection class, which provides functionality to:

- Crop model data to a target city or region.
- Generate urban masks based on the model urban and land-sea fraction.
- Generate rural masks based on land-cover and elevation criteria.
- Visualize the resulting masks.


### **urmask/utils.py**

Utility functions supporting data preprocessing, spatial operations, and visualization.


### **urmask/environment.txt**

Lists the main Python package dependencies and version requirements needed to reproduce the computational environment used by this repository.

Install all dependencies using:

pip install -r environment.txt


## Urban Boundary Datasets

The package can be used with different administrative boundary datasets, including:

### **Urban Audit dataset (URAU)**

Provides city (C) and functional urban area (F) polygons as defined by the EC-OECD city definition. This dataset is used for the Eurostat Urban Audit data collection European Commission, Eurostat, 2025). Used here as it is suitable for climate-service applications that require policy-relevant spatial definitions at the EU, national, and local government scales.

The data is distributed via https://gisco-services.ec.europa.eu/distribution/v2/urau/

_urban_rural_mask_example.ipynb_ uses 

**ref-urau-2024-100k.shp/**

**└── URAU_RG_100K_2024_4326.shp** 

downloaded from https://gisco-services.ec.europa.eu/distribution/v2/urau/urau-2024-files.html

### **Global Human Settlement - Urban Centre Database (GHS-UCDB)**

Provides globally consistent "urban centre" polygons developed by the European Commission, supporting global monitoring of policy frameworks and providing data for urban studies. GHS-UCDB inputs data from the GHS Layer (GHSL) and follows the “degree of urbanisation” (DEGURBA) methodology (Mari Rivero et al., 2026).

The data is distributed via https://human-settlement.emergency.copernicus.eu/ghs_stat_ucdb2015mt_r2019a.php

_urban_rural_mask_example.ipynb_ uses 

**GHS_STAT_UCDB2015MT_GLOBE_R2019A/**

**└── GHS_STAT_UCDB2015MT_GLOBE_R2019A_V1_2.gpkg**


downloaded from https://human-settlement.emergency.copernicus.eu/ghs_stat_ucdb2015mt_r2019a.php


## Citation

If you use urmask, please cite:

Puri, R., Teichmann, C., & Rechid, D. (2026). urmask: Urban and Rural Area Classification for High-Resolution Climate Model Data. Zenodo. https://doi.org/10.5281/zenodo.21341935


## References

- Diez-Sierra, J., Quintana, Y., Langendijk, G.S. et al. (2025). A global CORDEX-based dataset delineating urban areas and their surroundings to assess climate change in megacities. Scientific Data, 12, 1961. https://doi.org/10.1038/s41597-025-06257-1
  
- Diez-Sierra, J., Quintana, Y., Langendijk, G.S. et al. (2025). URCLIMASK: A Python Package for Delineating Urban Areas and Their Surrounding Reference Rural Regions from Regional Climate Models (RCMs) (Version v1.1.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.17257445
  
- European Commission. Eurostat. (2025). Urban Audit database [Dataset]. https://ec.europa.eu/eurostat/web/regions-and-cities

- Marí Rivero, I., Melchiorri, M., Florio, P. et al. (2026). GHS-UCDB R2024A: GHS Urban Centre Database 2025 [Dataset]. European Commission, Joint Research Centre. http://data.europa.eu/89h/1a338be6-7eaf-480c-9664-3a8ade88cbcd

- Sieck, K., Pinto, J. G., Geyer, B. et al. (2026). NUKLEUS – A First Kilometre Scale Multi-model Climate Ensemble for Germany: Evaluation. EGUsphere [preprint]. https://doi.org/10.5194/egusphere-2026-1024

