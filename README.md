# Urban Area Classification for Climate Model Data

This repository contains the **urmask** Python package, designed to derive urban and, optionally, rural masks from climate model data. Building on the **urclimask** Python package (Diez-Sierra et al., 2025), this **urmask** package extends the methodological framework for identifying urban and rural areas to convection-permitting regional climate models (CPRCMs).

**urmask** leverages each model's individual land representation to define urban and rural areas. This ensures that the resulting masks are consistent with the model's underlying surface characteristics.

The methodology uses the following model variables:

- Urban fraction – proportion of urbanized area within each grid cell.
- Land-sea mask – percentage of land area within each grid cell.
- Orography – surface elevation used to improve rural area selection.

To support climate-service applications, the methodology is additionally constrained by administrative boundaries, which can be adapted to different use cases. For example, users may define urban areas using City (C) or Functional Urban Areas (FUA) boundaries from the Urban Audit (URAU) dataset. This approach enables the generation of administratively meaningful urban and rural masks that remain consistent with the land-cover representation of individual climate models.

## Repository Contents

urban_rural_mask_example.ipynb

A Jupyter Notebook demonstrating the complete workflow for generating urban and rural masks using:

A regional climate model land-representation dataset.
An administrative city boundary polygon.

Administrative boundary datasets are not distributed with this repository and must be downloaded separately from the sources listed below.

urmask/
urmask.py

Contains the UrbanRuralSelection class, which provides functionality to:

Crop model data to a target city or region.
Generate urban masks based on model urban fraction.
Generate rural masks based on land-cover and elevation criteria.
Visualize the resulting masks.
utils.py

Utility functions supporting data preprocessing, spatial operations, and visualization.

environment.txt

Lists the main Python package dependencies and version requirements needed to reproduce the computational environment used by this repository.

Install all dependencies using:

pip install -r environment.txt

Usage

See urban_rural_mask_example.ipynb for a step-by-step demonstration of the urban and rural mask generation workflow.

Administrative Boundary Datasets

The package can be used with different administrative boundary datasets, including:

Urban Audit (URAU)

Provides harmonized city and Functional Urban Area boundaries across Europe, suitable for climate-service applications requiring policy-relevant spatial definitions.

Urban Centre Database (UCDB)

Provides globally consistent urban centre boundaries and associated metadata, enabling application of the methodology beyond Europe.

Purpose

The primary objective of urmask is to provide a reproducible and flexible framework for defining urban and rural areas in high-resolution climate simulations that:

Preserves consistency with the climate model's own land-surface representation.
Produces administratively meaningful spatial masks.
Supports climate-impact assessments and climate-service development.
Enables comparable analyses across different regional climate models and urban boundary definitions.



**- urban_rural_mask_example.ipynb**
_A Jupyter notebook demonstrating how to calculate urban and rural areas using a land representation dataset file for a regional climate model and administrative city polygon (not included, please download from the sources below)._

- **urmask**
    - **_urmask.py**
        _A Python module that includes the UrbanRuralSelection class that crops the city area, defines urban and rural masks, and gives a plotting function for the same._

    - **_utils.py**

      Note: shapefiles for URAU and UCDB are available under:
        

- **environment.txt**
  _A text file listing the main Python packages required and their versions to reproduce the computational environment needed to run the two files above. You can install these using pip install -r environment.txt._

## Usage
See urban_rural_mask_example.ipynb for a step-by-step example.

### Shapefiles

- **Urban Audit Dataset**
- **Urban Center Database**


