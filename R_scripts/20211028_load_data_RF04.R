nc_path = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/Outputs/"
rr.stack <- stack(file.path(nc_path, "LES_outputs_rasters/RF04/rr.stack.nc"))
wind.stack <- stack(file.path(nc_path, "LES_outputs_rasters/RF04/wind.stack.nc"))
psfc.stack <- stack(file.path(nc_path, "LES_outputs_rasters/RF04/psfc.stack.nc")) * 0.01 # convert from Pa to hPa 

current.date = "2021-07-30"

nc_file = paste(nc_path, "LES_outputs_rasters/RF04/RF04_new.nc", sep = "/")
nc_data <- nc_open(nc_file)
nc_times <- ncvar_get(nc_data, "Times")
Times <- as.POSIXct(strptime(nc_times, "%Y-%m-%d_%H:%M:%S"), tz = "UTC")