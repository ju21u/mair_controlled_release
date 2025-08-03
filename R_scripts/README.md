# R_scripts Folder

> ⚠️ **Note:** Many scripts refer to specific file names or hardcoded paths at the beginning.  
> You will need to update these paths and file references to match your own directory structure and data locations.

This directory contains R scripts used for data loading, pre-processing, and methane emission estimation for the MethaneAIR Controlled Release campaign.

---

## Data Loading & Cleaning

| Script                          | Purpose                                                                 |
|---------------------------------|-------------------------------------------------------------------------|
| 20211028_load_data_RF04.R       | Loads data for research flight RF04. Likely includes. |
| 20211101_load_data_temp.R       | Load temporary; serve as a staging area. |
| 20220111_clean_extents.R        | Cleans or adjusts spatial extents for analysis (e.g., bounding boxes or rasters). |
| 20230207_clean_1x1_flare_release.R | Output estimated emissions of the controlled releases and the flares (to the left of the controlled release site in RF04 and RF05). |
| 20230328_clean_1x1_releaseE.R   | Output estimated emissions of the RF01E and RF03E controlled releases. |

---

## Emission Estimation

| Script                    | Purpose                                                                 |
|---------------------------|-------------------------------------------------------------------------|
| estimate_emissions_v11.R | Most recent version of the main emissions estimation function. Uses cropped scenes, wind, and pressure data. |
| estimate_emissions_v8.R  | Older version of the above; retained for historical comparison. Superseded by `v11`. |

---

## Shared Resources

| Script                    | Purpose                                                                 |
|---------------------------|-------------------------------------------------------------------------|
| load_functions_constants.R | Defines custom functions and physical or campaign-specific constants. Used across scripts. |
| load_libraries.R         | Loads all required R packages. Run this script before others. |
| load_required_files.R    | Loads common external input files such as raster templates, masks, or meteorological inputs. |
