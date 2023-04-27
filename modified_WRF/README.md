# Level 4 Large Eddy Simulation (LES) using [WRF](https://www.mmm.ucar.edu/weather-research-and-forecasting-model)

## Data Version Control (DVC)

We will be using DVC for version controlling large files. Please see the following [readme](https://github.com/methanesat-org/methanesat/blob/main/DVC.md) for setup instructions. 

Also the [dvc.org](dvc.org) website has very well written documentation

## Modified WRF

This repo contains the modified files for the WRF-Chem version 3.9.1.1. The code was originally modified by Lauvaux et al., 2012 [link here](https://acp.copernicus.org/articles/12/337/2012/). 

### Docker

A docker container has been created to encapsulate building/running WRF v3.9.1 with these modified files. The dockerfile can be found in the [msat-platform repo](https://github.com/methanesat-org/msat-platform/tree/main/build/docker)

### Manually

To use this code manually, first download the 3.9.1.1 model from [here (.zip)](https://github.com/NCAR/WRFV3/archive/V3.9.1.1.zip) or [here (.tar.gz)](https://github.com/NCAR/WRFV3/archive/V3.9.1.1.tar.gz).

Then integrate these modified files to the model that you just downloaded. 

Besides the WRF folder, you also need the WRF Preprocessing System, which can be downloaded [here](https://github.com/NCAR/WPS)
Note that you still need the met files (HRRR from [here](http://home.chpc.utah.edu/~u0553130/Brian_Blaylock/cgi-bin/hrrr_download.cgi?model=hrrr&field=prs&date=2020-04-05&link2=grib2)) and WPS_GEOG ([here](https://www2.mmm.ucar.edu/wrf/src/wps_files/geog_complete.tar.gz) and [here](https://www2.mmm.ucar.edu/wrf/src/wps_files/geog_new3.9.tar.bz2))


You may need to downloaded a Vtable for HRRR met files if you want to use HRRR as your driven meterological data (can be downloaded [here](https://github.com/blaylockbk/Ute-WRF-User-Group/blob/master/Blaylock/WPS_for_HRRR/Vtable.HRRR.bkb)). 




##### Additional Instructions for ERA5 data

*** Note that the files need to be subsetted in the operational scale ****

*** For now, we download all the files for the entire globe ***


- Download surface and pressure data from ucar using the load_ERA5.sh template in WPS folder
- Use the Vtable.ECMWF to ungrib these data. To do so, @ WPS folder, ln -sf ungrib/Variable_Tables/Vtable.ECMWF Vtable before running ungrib.exe 
- Make sure that the resolutions are correct. The ERA5 namelist on this branch should be correct. 
- Follow the remaining steps as usual 
