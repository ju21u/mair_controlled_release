# Load data MethaneAIR Data 
#base.dir = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/"
#source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/Bella/20211001_load_libraries.R")
#load("/Users/apisada/Documents/R_scripts/heat.275.RData")


nc_path = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/Outputs/LES_outputs_rasters/RF05"
rr.stack <- stack(file.path(nc_path, "rr.stack.nc"))
wind.stack <- stack(file.path(nc_path, "wind.stack.nc"))
psfc.stack <- stack(file.path(nc_path, "psfc.stack.nc")) * 0.01 # convert from Pa to hPa
load(file.path(nc_path, 'saved_variables.RData'))


#nc_path = "/Users/apisada/Google\ Drive/My\ Drive/Research/Harvard_Research/MethaneAIR"

#summary.steve <- read.csv(file.path(nc_path, "From_Steve/all.segs/RF05/5x1/RELEASE_RASTERS/RF05_release_summary.csv"), sep = " ")

