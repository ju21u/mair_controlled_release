# go to the SEGMENT_RASTER folder
# FIX THE PRESSURE --> no pressure from Steve atm 
cat('\014')
graphics.off()
rm(list = ls())
old.dir = getwd()

source("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/estimate_emissions_v11.R")
source("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/unique_scripts/20220111_clean_extents.R")
source("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/common_scripts//load_libraries.R")
source("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/common_scripts/load_functions_constants.R")
source("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/R_scripts/MethaneAIR/common_scripts/load_required_files.R")

export_data = FALSE
resolution = "1x1"
crop.radius = 700 # 0.75 km 
N_i = 1000 # number of iterations
out.folder = file.path("/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/Outputs/" ) #CSVs", resolution, paste(as.character(crop.radius), "m", sep = ""))
cropped_plume_folder = "/n/holylfs04/LABS/wofsy_lab/Users/achulakadabba/MethaneAIR/controlled_plume_rasters"
########################################
########################################
########################################

# current.research.flight = "RF05"
# source.name = "flare"
########################################
########################################
########################################

raw_bootstrapped_dataframe = data.frame(RF = NULL, seg = NULL, formatted.time = NULL, rawbt = NULL)
  
MASS_WINDS_release = NULL
MASS_WINDS_flare = NULL

bigbox = extent( c( -102.3172,-102.26299, 32.03887, 32.10841 ) )

#thres_all  = 2 ## number of sd's above the median
thres_all = c(1.5) #, 2) # , 2.5, 3) #c(1, 1.5, 2, 2.5, 3)
for (current_thres in thres_all){
  min_pixels_in_a_clump = 20 ## number of points in a clump (> 2000 sq m ++) 
  
  output.pdf = paste(Sys.Date(), "bg_summary", resolution, as.character(crop.radius/100),"km.thres", current_thres, "n.pixels", min_pixels_in_a_clump, '.pdf', sep = '.')
  # 'summary_entire_plume_10m_10m_1km_both'
  # DEFINE EXTENTS FOR THE SOURCES OF INTEREST 
  e = extent( c( -102.30872,-102.27299, 32.04887, 32.09841 ) )
  
  
  RF04_flare_summary<- data.frame(matrix(vector(), 0, 16,
                                         dimnames=list(c(),  c("SEGMENT", "IME EMISS", "IME AREA", "IME IME", "clumps", "RATIO EMISS", "RATIO SD", "RATIO 25",  "RATIO 975", "IME 25", "IME 975", "SD","ns NOISE", "bg NOISE"," WIND", "THRES") )),
                                  stringsAsFactors=F)
  RF04_release_summary <- RF04_flare_summary
  RF05_flare_summary <- RF04_flare_summary
  RF05_release_summary <- RF04_flare_summary
  
  RF04_release_mi_summary <- data.frame(matrix(vector(), 0, 16,
                                               dimnames=list(c(),  c("center.of.mass.lon", "center.of.mass.lat", "angle.deg", "mass", "area", "plume.length.scale", "moments.minor", "moments.major", "eccentricity", "n.clumps.in.plume", "sig.sel", "min.px.in.clump", "inflow", "inflow.sd", "threshold", "seg") )), #"n.clumps.in.plume", "sig.sel", "inflow",  "RATIO 975", "Kg") )),
                                        stringsAsFactors=F)
  RF05_release_mi_summary <- RF04_release_mi_summary
  RF04_flare_mi_summary <- RF04_release_mi_summary
  RF05_flare_mi_summary <- RF04_release_mi_summary
  

  # # RF04
  ###### COMMENTS #################
  segment.inputs = seq(9)#@c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21) #seq(9) #(1, 5) #, 10, 12, 14)


  static.release.inflow = c(-102.3060, -102.300, 32.0503, 32.052)
  static.release = c(-102.3045, -102.3005, 32.052, 32.06000) #c(-102.3045, -102.298, 32.052, 32.06000) #
  static.release.ext <-  data.frame(t(static.release))
  static.release.plot.lon <- static.release.inflow[1:2]
  static.release.plot.lat <- static.release.inflow[3:4]

  colnames(static.release.ext) <- c("xmin", "xmax", "ymin", "ymax")


  static.flare.inflow = c(-102.3045, -102.300, 32.0503, 32.052)
  static.flare = c(-102.294, -102.284, 32.053, 32.075)
  static.flare.ext <-  data.frame(t(static.flare))
  colnames(static.flare.ext) <- c("xmin", "xmax", "ymin", "ymax")
  static.flare.plot.lon <- static.flare.inflow[1:2]
  static.flare.plot.lat <- static.flare.inflow[3:4]
  RF04_n_pixel_summary = c()
  RF05_n_pixel_summary = c()
  RF04_all_conc_summary = list()
  RF05_all_conc_summary = list()
  RF04_small_conc_summary = list()
  RF05_small_conc_summary = list()
  RF04_large_conc_summary = list()
  RF05_large_conc_summary = list()
  RF05_inflow_summary = list()
  RF04_inflow_summary = list()
  #
 # for (current.idx in c(7)) {
  for (current.idx in seq(length(segment.inputs))) {
  #  if(!is.na(sum(RF04.flare.inflows[current.idx, ]))) {
    current.research.flight  = "RF04"
    #current.path = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF04/1x1/SEGMENT_RASTERS/"
    current.segment.input = segment.inputs[current.idx]

    # load what needs to be loaded for the estimate_emissions function to work
      load_files_outputs = load_files(current.research.flight = current.research.flight, res = resolution)
      source(load_files_outputs[[1]])
      current.path = load_files_outputs[[2]]
      setwd(current.path)

      out.all = estimate_emissions( current.research.flight = current.research.flight,
                      #xy.input = xy.flare, box.inflow.input = RF04.release.inflows[current.idx, ],
                      xy.input = xy.flare, box.inflow.input = static.flare.inflow,
                      e.input = extent(static.flare), source.name.input = 'flare',
                      #e.input = extent(RF04.release[current.idx, ]), source.name.input = 'release',
                      current.segment.input  = current.segment.input, res = resolution, N_i = N_i,
                      plume_crop = NULL, save_raw = TRUE)
      out.row = out.all[[1]]
      out_n_pixels = out.all[3]
      RF04_n_pixel_summary <- rbind(RF04_n_pixel_summary, out_n_pixels)
      RF04_flare_summary <- rbind(RF04_flare_summary, out.row)
      RF04_all_conc_summary <- append(RF04_all_conc_summary, out.all[4])
      RF04_small_conc_summary <- append(RF04_small_conc_summary, out.all[5])
      RF04_large_conc_summary <- append(RF04_large_conc_summary, out.all[6])
      RF04_inflow_summary <- append(RF04_inflow_summary, out.all[8])
     # browser()

     # write.csv(out.all[[2]][[4]], file.path(out.folder, paste(format(Sys.Date(), "%Y%m%d"),'_', current.research.flight, '_', current.segment.input, '.csv', sep = '' )))

  }
  write.csv(RF04_flare_summary,file.path(out.folder, paste(format(Sys.Date(), "%Y%m%d"),'_RF04_flare_1000_',  as.character(current_thres), '.csv', sep = '' )))

  e_current_background = extent(c(-102.32, -102.31, 32.10, 32.11))

  ###### COMMENTS #################
  #ser()
#   #RF05
  ########################################
  segment.inputs =  c(1, 3, 5, 7, 9, 11,  13, 15, 17, 19, 21, 23, 25)  #c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22)  #c(1, 3, 5,  7, 9, 11,  13, 15, 17, 21, 23)  #c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23) #, 25, 27, 29) #c( 25, 27, 29) # c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23,
  
  
  static.flare.inflow = c(-102.3045, -102.300, 32.046, 32.049)
  static.flare = c(-102.294, -102.284, 32.0505, 32.075)
  
for (current.idx in seq(length(segment.inputs))) {
    current.research.flight  = "RF05"
    current.segment.input = segment.inputs[current.idx]
    
    load_files_outputs = load_files(current.research.flight = current.research.flight, res = resolution)
    source(load_files_outputs[[1]])
    current.path = load_files_outputs[[2]]
    setwd(current.path)

    if(!is.na(sum(RF05.release.inflows[current.idx, ]))) {

     out.all = estimate_emissions( current.research.flight = current.research.flight,
                                    xy.input = xy.flare, box.inflow.input = static.flare.inflow,
                                   e.input = extent(static.flare), source.name.input = 'flare',
                                   # e.input = extent(RF05.release[current.idx, ]), source.name.input = 'release',
                                    current.segment.input  = current.segment.input, 
                                   res = resolution, N_i = N_i, save_raw = TRUE)
      out.row = out.all[[1]]
      out_n_pixels = out.all[3]
      RF05_n_pixel_summary <- rbind(RF05_n_pixel_summary, out_n_pixels)
      RF05_all_conc_summary <- rbind(RF05_all_conc_summary, out.all[4])
      RF05_small_conc_summary <- rbind(RF05_small_conc_summary, out.all[5])
      RF05_large_conc_summary <- append(RF05_large_conc_summary, out.all[6])
      RF05_inflow_summary <- append(RF05_inflow_summary, out.all[8])
      MASS_WINDS_release = rbind(MASS_WINDS_release, out.all[[2]])
      RF05_flare_summary <- rbind(RF05_flare_summary, out.row)}
    print(paste("FINISHED SEG", as.character(segment.inputs[current.idx])))
   # write.csv(out.all[[2]][[4]], file.path(out.folder, paste(format(Sys.Date(), "%Y%m%d"),'_', current.research.flight, '_', current.segment.input, '.csv', sep = '' )))
    
  }
  write.csv(RF05_flare_summary,file.path(out.folder, paste(format(Sys.Date(), "%Y%m%d"), '_RF05_flare_1000_',  as.character(current_thres), '.csv', sep = '' )))
  
  print("=====================================")
  print("=====================================")
  print(paste(" DONE WITH CURRENT THRES:", as.character(current_thres)))
  print("=====================================")
  print("=====================================")

}
