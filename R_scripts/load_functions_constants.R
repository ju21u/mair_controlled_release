# Functions and constants 

## Constants ## ---------------------------------------------------
xy.flare = list("x" =  -102.2875 , "y" = 32.056)
xy.release = list("x" =  -102.300687 , "y" =  32.053111)

Earth_radius = 6378137 # meters

# load color palette 
script_folder = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/common_scripts/"
load(file.path(script_folder, "heat.275.RData"))

source(file.path(script_folder, "functions/moments.fcn2.r"))
## functions ## ---------------------------------------------------
molef.to.kg = function( XC, Psurf, mole.wt=.016) {
  # H = 1.38e-23*Tsurf/( mole.wt/6.03e23 * 9.8)
  N.CH4 = Psurf*100 / 9.8  * XC * mole.wt/ 0.029
  #       hPa  - Pa   g     CMDMF
  return( N.CH4)
} 

generate_binary_array = function(input_array, conditions){
  binary_array <- array(0, c(dim(input_array)))   
  binary_array[conditions] = 1
  
  return(binary_array)
}
find_wings = function (new_slope, x, y, wing_size){
  theta_rad = atan2(x*new_slope, x)
  x_delta = cos(theta_rad) * wing_size
  y_delta = sin(theta_rad) * wing_size
  x_end_1 = x - x_delta
  x_end_2 = x + x_delta
  y_end_1 = y - y_delta
  y_end_2 = y + y_delta
  out_array = list(list(x_end_1, x_end_2), list(y_end_1, y_end_2))
  return(out_array) 
  }

rms.cal = function(input){
  output = sqrt(mean(input**2))
  return(output)
}

matrix_rotation = function(x, y, z, theta){
  x_new = x*cos(theta) - y*sin(theta) 
y_new = x*sin(theta) + y*cos(theta) 
x_new_positive = x_new + abs(min(x_new)) 
y_new_positive = y_new + abs(min(y_new))
return(list(x_new_positive, y_new_positive, z))}

find_y = function(xy_fit, x){
  y = xy_fit[1] * x + xy_fit[2]
  return (y)
}

perpen_cross = function(slope, interp){
  new_slope = -1/slope
  return(new_slope)
}
         
IME.cal = function(t_str = 349, t_end = 359, x.range = 58:62, y.range = 58:62, pass_section = "6_12", source.type = "flr") {
  idx = t_str:t_end
  #IME = total_observed_mass * area
  current_summary = paste("summary.2_", pass_section, sep = "")
  wind.2 = U_surface[58:62, 58:62, idx]**2 +  V_surface[58:62, 58:62, idx]**2
  if(length(dim(wind.2)) > 2){
    wind = sqrt(colSums(wind.2, dims = 2))}
  else{wind = sqrt(sum(wind.2))}
  mass_str = paste("total.mass.", source.type, ".2", sep = "")
  IME =  Summary_J.2[ current_summary, mass_str] 
  U_eff = log(wind) + 0.6 #m/s
  area_str =  paste("area.", source.type, ".m2.2", sep = "")
  L = sqrt(Summary_J.2[current_summary, area_str])
  EMISS = IME *  U_eff / L #kg/s
  EMISS_hr = EMISS * 3600
  if(is.nan(EMISS_hr)){EMISS_hr = 0}
  return( EMISS_hr)}
IME.cal.simple = function(t_str = 349, t_end = 359, x.range = 58:62, y.range = 58:62, ime_input = 0, area_input = 0) {
  idx = t_str:t_end
  #IME = total_observed_mass * area
  #current_summary = paste("summary.2_", pass_section, sep = "")
  wind.2 = U_surface[x.range, y.range, idx]**2 +  V_surface[x.range, y.range, idx]**2
  if(length(dim(wind.2)) > 2){
    wind = sqrt(colSums(wind.2, dims = 2))}
  else{wind = sqrt(sum(wind.2))}
  #mass_str = paste("total.mass.", source.type, ".2", sep = "")
  IME =  ime_input #Summary_J.2[ current_summary, mass_str] 
  U_eff = log(wind) + 0.6 #m/s
  #area_str =  paste("area.", source.type, ".m2.2", sep = "")
  L = sqrt(area_input)
  EMISS = IME *  U_eff / L #kg/s
  EMISS_hr = EMISS * 3600
  if(is.nan(EMISS_hr)){EMISS_hr = 0}
  return( EMISS_hr)}

emission.cal.IME = function(wind = 0,ime_input = 0, area_input = 0, coef = 1) {
  #mass_str = paste("total.mass.", source.type, ".2", sep = "")
  IME =  ime_input #Summary_J.2[ current_summary, mass_str] 
  U_eff = coef*log(wind) + 0.6 #m/s
  #area_str =  paste("area.", source.type, ".m2.2", sep = "")
  L = sqrt(area_input)
  EMISS = IME *  U_eff / L #kg/s
  EMISS_hr = EMISS * 3600
  if(!is.numeric(EMISS_hr)){EMISS_hr = 0}
  return( EMISS_hr)}

emission.cal.coef.IME = function(wind = 0, ime_input = 0,  area_input = 0, les_ime_input = 0, les_area_input = 0, les_emiss_hr = 988) {
  L = sqrt(area_input) # scalelength from the observation 
  IME = ime_input # IME from the observation 
  les.emiss = les_emiss_hr/3600
  les.L = sqrt(les_area_input)
  new.U.eff = les.emiss / les_ime_input * les.L # [kg/s]/[kg] * [m] == [m/s]
  EMISS = IME *  new.U.eff / L #kg/s
  EMISS_hr = EMISS * 3600
  if(is.nan(EMISS_hr)){EMISS_hr = 0}
  # print(paste("predicted LES emiss", as.character(les_ime_input * new.U.eff / les.L * 3600)))
  # print(paste("new u_eff: ", as.character(new.U.eff)))
  # print(paste("new les.emiss: ", as.character(les.emiss)))
  # print(paste("new les_ime_input: ", as.character(les_ime_input)))
  # print(paste("new les.L: ", as.character(les.L)))
  return( EMISS_hr)}


get.dir = function(U = 3, V = 3){
  slope = mean(V/U)
  rads = atan(slope)
  degs = rads * 180/pi
  return(list(slope, degs, rads)) # in degrees and radians 
}
# select.area = function(lat, lon, R, ws.dir){
#   
#   return(out.matrix)
# }


plot_xyz = function(x, y, z, title){
  Col <- heat.275[as.numeric(cut(z, breaks = 275))]
  plot(x, y, pch = 16,  col = Col, main = title)
  image.plot(legend.only=T, zlim=range(z), col = heat.275)
}


create.raster.stack = function (dframe, dlong, dlat, dtime, outfile){
  
  rr.stack <- raster(xmn= min(dlong), xmx = max(dlong),
                     ymn= min(dlat), ymx = max(dlat),
                     crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
  
  for(t_idx in seq(dim(dtime))){
    rrr = dframe[, , 1 ,t_idx] #surface wind only
    rr.0 = t( rrr) [ ncol( rrr ):1, ]
    Rrr = raster( rr.0, xmn= min(dlong), xmx = max(dlong),
                  ymn= min(dlat), ymx = max(dlat),
                  crs = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" )
    rr.stack <- stack(rr.stack, Rrr)
  }
  names(rr.stack) <-seq(nlayers(rr.stack))
  dt <- as.data.frame(dtime) 
  outvar <- setZ(rr.stack, dt, "Times")
  #outvar <- setZ(rr.stack, dt, "Times")
  raster::writeRaster(outvar, filename = outfile, format="CDF", overwrite=TRUE)
  return(outvar)
}


ueff.cal = function(les_ime_input = 0, les_area_input = 0, les_emiss_hr = 988) {
  # IME from the LES was given
  les.emiss = les_emiss_hr/3600
  les.L = sqrt(les_area_input)  # scalelength from LES
  new.U.eff = les.emiss / les_ime_input * les.L
  return(new.U.eff)
}


