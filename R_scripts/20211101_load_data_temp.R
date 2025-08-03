# Load data MethaneAIR Data 
nc_path = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/Outputs/LES_outputs_rasters/RF05"
rr.stack <- stack(file.path(nc_path, "rr.stack.nc"))
wind.stack <- stack(file.path(nc_path, "wind.stack.nc"))
psfc.stack <- stack(file.path(nc_path, "psfc.stack.nc")) * 0.01 # convert from Pa to hPa
load(file.path(nc_path, 'saved_variables.RData'))
