# MethaneAIR Controlled Release Experiment 
MethaneAIR's Blinded Volume Controlled Release Experiments in 2021

# MethaneAIR Controlled Release Study (2025 Updates)

This repository contains scripts and configurations for the MethaneAIR Controlled Release campaign, including emission quantification, data cleaning, and WRF-STILT model setup.

##  Repository Structure

```
mair_controlled_release-2025_updates/

 R_scripts/                  # Main R scripts for data processing and IME emission estimation
 DI_scripts/                 # Scripts for Divergent Integral methods
 modified_WRF/               # Customized WRF & WPS setup
    WRF/
    WPS/
 LICENSE
 README.md                  # You're reading it
 .gitignore
```

##  Getting Started

### Prerequisites

- R (>= 4.0)
- Required R packages (install via `load_libraries.R`)
- WRF model build environment (e.g., Intel compiler, netCDF, MPI)

### Setup

```bash
git clone https://github.com/ju21u/mair_controlled_release-main.git
cd mair_controlled_release-main

# In R:
source('R_scripts/load_libraries.R')
source('R_scripts/load_required_files.R')
```

### Running Emission Estimation

```r
source('R_scripts/estimate_emissions_v11.R')
```

##  Data

Large `.nc` files and model output are not tracked in Git. Please place external data in a `data/` directory as needed.

##  WRF Setup

The `modified_WRF/` folder contains a customized WRF-WPS setup for processing ERA5/GDAS data. See the README inside that folder for WPS usage.

##  Tips for Maintenance

- Use consistent naming (avoid timestamps in filenames).
- Archive outdated versions instead of keeping many copies.
- Prefer `load_libraries.R` to manage package dependencies.

##  License

This project is licensed under the terms of the LICENSE file in the root directory.

## Visual Output

Example plume detection result:

![Plume Visualization](https://user-images.githubusercontent.com/31904333/234888223-c1c171b5-d9c8-426f-bc72-3c479c5d2e02.png)

## Data and Resources

Additional files, including large datasets and processed outputs, are available in the shared Google Drive:

[MethaneAIR Shared Drive](https://drive.google.com/drive/folders/1Xg57yA2dFUpWyzNnZoe-Px13pd4NpnWq?usp=share_link)
