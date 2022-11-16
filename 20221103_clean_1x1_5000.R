# go to the SEGMENT_RASTER folder
# FIX THE PRESSURE --> no pressure from Steve atm 
cat('\014')
graphics.off()
rm(list = ls())
old.dir = getwd()
library(dplyr)
library(bootstrap)
#source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF04/load_libraries.R")
#source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/Bella/20211001_load_libraries.R")
#ource("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean/estimate_emissions.R")
#source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean/20220111_clean_extents.R")
source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean/estimate_background.R")
source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean/load_libraries.R")

source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean/20220111_clean_extents.R")
source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean/load_functions_constants.R")



resolution = "1x1"
crop.radius = 5000 # 5km 

out.folder = file.path("/Volumes/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/Outputs/CSVs", resolution, paste(as.character(crop.radius), "m", sep = ""))
########################################
########################################
########################################

# current.research.flight = "RF05"
# source.name = "flare"
########################################
########################################
########################################
load("/Users/apisada/Documents/R_scripts/heat.275.RData")

MASS_WINDS_release = NULL
MASS_WINDS_flare = NULL

#en0  = 2 ## number of sd's above the median
en0_all = c(1.5) #c(1, 1.5, 2, 2.5, 3) #c(1, 1.5, 2, 2.5, 3)
for (en0 in en0_all){
  n.sel = 20 ## number of points in a clump (> 2000 sq m ++) 
  
  output.pdf = paste(Sys.Date(), "bg_summary", resolution, as.character(crop.radius/100),"km.thres", en0, "n.pixels", n.sel, '.pdf', sep = '.')
  # 'summary_entire_plume_10m_10m_1km_both'
  # DEFINE EXTENTS FOR THE SOURCES OF INTEREST 
  e = extent( c( -102.30872,-102.27299, 32.04887, 32.09841 ) )
  #e.rel = extent ( c( -102.3086, -102.2981, 32.05030, 32.06312))
  
  
  
  
  ##  define the inflow box
  
  
  #################
  #################
  
  #9999 #300 #1000
  # pdf(file = file.path("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Outputs/PDFs/", output.pdf),
  #     width = 7, # The width of the plot in inches
  #     height = 8) # The height of the plot in inches

  # graphics.off()
  
  xy.flare = list("x" =  -102.2875 , "y" = 32.056)
  xy.release = list("x" =  -102.300687 , "y" =  32.053111)
  
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
  
  
  # output.pdf = paste(Sys.Date(), 'summary', en0, n.sel, 'pdf', sep = '.')
  # pdf(file = file.path("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Outputs/PDFs/", output.pdf),
  #     width = 7, # The width of the plot in inches
  #     height = 8) # The height of the plot in inches
  
  
  # RF04
  segment.inputs = seq(9) #(1, 5) #, 10, 12, 14)
  
  static.release.inflow = c(-102.3060, -102.300, 32.0503, 32.052)
  static.release = c(-102.3045, -102.3005, 32.052, 32.06000)
  static.release.ext <-  data.frame(t(static.release)) 
  static.release.plot.lon <- static.release.inflow[1:2]
  static.release.plot.lat <- static.release.inflow[3:4]
  
  colnames(static.release.ext) <- c("xmin", "xmax", "ymin", "ymax")
  
  
  static.flare.inflow = c(-102.3045, -102.300, 32.0503, 32.052)
  static.flare = c(-102.294, -102.284, 32.055, 32.075)
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
  for (current.idx in seq(length(segment.inputs))) {
  #  if(!is.na(sum(RF04.flare.inflows[current.idx, ]))) {
    current.research.flight  = "RF04"
    current.path = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF04/1x1/SEGMENT_RASTERS/"
    current.segment.input = segment.inputs[current.idx]

     # out.all = estimate_emissions( current.research.flight = current.research.flight ,
     #                  xy.input = xy.flare, box.inflow.input = RF04.flare.inflows[current.idx, ],
     #                  #xy.input = xy.flare, box.inflow.input = static.flare.inflow,
     #                  #e.input = extent(static.flare), source.name.input = 'flare',
     #                  e.input = extent(RF04.flare[current.idx, ]), source.name.input = 'flare',
     #                  current.segment.input = current.segment.input, res = resolution)
     #  out.row = out.all[[1]]
     #  RF04_flare_summary <- rbind(RF04_flare_summary, out.row)
    #  }
      setwd(current.path)
     

   # if(!is.na(sum(RF04.release.inflows[current.idx, ]))) {

      out.all = estimate_emissions( current.research.flight = current.research.flight,
                      xy.input = xy.release, box.inflow.input = RF04.release.inflows[current.idx, ],
                      #xy.input = xy.release, box.inflow.input = static.release.inflow,
                      e.input = extent(static.release), source.name.input = 'release',
                      #e.input = extent(RF04.release[current.idx, ]), source.name.input = 'release',
                      current.segment.input  = current.segment.input, res = resolution)
      out.row = out.all[[1]]
      out_n_pixels = out.all[3]
      RF04_n_pixel_summary <- rbind(RF04_n_pixel_summary, out_n_pixels)
      RF04_release_summary <- rbind(RF04_release_summary, out.row)
      RF04_all_conc_summary <- append(RF04_all_conc_summary, out.all[4])
      RF04_small_conc_summary <- append(RF04_small_conc_summary, out.all[5])
      RF04_large_conc_summary <- append(RF04_large_conc_summary, out.all[6])
      RF04_inflow_summary <- append(RF04_inflow_summary, out.all[8])
      #browser()

  }
  
  e_current_background = extent(c(-102.32, -102.31, 32.10, 32.11))
  #browser()
   #}
  # write.csv(RF04_flare_summary,file.path(out.folder, paste('RF04_flare_summary_thres_', as.character(en0), '.csv', sep = '' )))
  # write.csv(RF04_release_summary,file.path(out.folder, paste('RF04_release_summary_thres_',  as.character(en0), '.csv', sep = '' )))
  # 
  # write.csv(RF04_flare_mi_summary,file.path(out.folder, paste('RF04_flare_mi_summary_thres_', as.character(en0), '.csv', sep = '' )))
  # write.csv(RF04_release_mi_summary,file.path(out.folder, paste('RF04_release_mi_summary_thres_',  as.character(en0), '.csv', sep = '' )))
  
  #ser()
  #RF05
  ########################################
  segment.inputs = c(1, 3, 5,  7, 9, 11,  13, 15, 17, 21, 23)  #c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23) #, 25, 27, 29) #c( 25, 27, 29) # c(1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23,
  # c(1, 11, 13, 15, 17, 19, 21, 23, 25, 3, 30, 5, 7, 9 )
  for (current.idx in seq(length(segment.inputs))) {
    #current.idx = current.idx  + 12 
    current.research.flight  = "RF05"
    current.path = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF05/1x1/SEGMENT_RASTERS/"
    current.segment.input = segment.inputs[current.idx]
    
    # print(as.character(segment.inputs[current.idx - 12]))
    if(!is.na(sum(RF05.flare.inflows[current.idx, ]))) {
      
      # out.all = estimate_emissions( current.research.flight = current.research.flight,
      #                               xy.input = xy.flare, box.inflow.input = RF05.flare.inflows[current.idx, ],
      #                               e.input = extent(RF05.flare[current.idx, ]), source.name.input = 'flare',
      #                               current.segment.input = current.segment.input, res = resolution)
      # 
      # out.row = out.all[[1]]
      # MASS_WINDS_flare = rbind(MASS_WINDS_flare, out.all[[2]])
      # print(segment.inputs[current.idx])
      # RF05_flare_summary <- rbind(RF05_flare_summary, out.row)
      # 
      # print("done flare")
    }
    
    
    if(!is.na(sum(RF05.release.inflows[current.idx, ]))) {
      
      out.all = estimate_emissions( current.research.flight = current.research.flight,
                                    xy.input = xy.release, box.inflow.input = RF05.release.inflows[current.idx, ],
                                    e.input = extent(RF05.release[current.idx, ]), source.name.input = 'release',
                                    current.segment.input  = current.segment.input, res = resolution)
      out.row = out.all[[1]]
      out_n_pixels = out.all[3]
      RF05_n_pixel_summary <- rbind(RF05_n_pixel_summary, out_n_pixels)
      RF05_all_conc_summary <- rbind(RF05_all_conc_summary, out.all[4])
      RF05_small_conc_summary <- rbind(RF05_small_conc_summary, out.all[5])
      RF05_large_conc_summary <- append(RF05_large_conc_summary, out.all[6])
      RF05_inflow_summary <- append(RF05_inflow_summary, out.all[8])
      MASS_WINDS_release = rbind(MASS_WINDS_release, out.all[[2]])
      RF05_release_summary <- rbind(RF05_release_summary, out.row)}
    print("done release")
    print(paste("FINISHED SEG", as.character(segment.inputs[current.idx])))
  
    # 
  }
  #RF05_release_summary
  #browser()
 
  # RF04 c(1, 5, 10, 12, 14)
  # RF05 c(1, 11, 13, 15, 17, 19, 21, 23, 25, 3, 30, 5, 7, 9 )
  # write.csv(RF05_flare_summary,file.path(out.folder, paste('RF05_flare_summary_thres_', as.character(en0),'.csv', sep = '' )))
  # write.csv(RF05_release_summary,file.path(out.folder, paste('RF05_release_summary_thres_', as.character(en0), '.csv', sep = '' )))
  # 
  # write.csv(RF05_flare_mi_summary,file.path(out.folder, paste('RF05_flare_mi_summary_thres_', as.character(en0), '.csv', sep = '' )))
  # write.csv(RF05_release_mi_summary,file.path(out.folder, paste('RF05_release_mi_summary_thres_',  as.character(en0), '.csv', sep = '' )))
  
  #browser()
  print("=====================================")
  print("=====================================")
  print(paste(" DONE WITH CURRENT THRES:", as.character(en0)))
  print("=====================================")
  print("=====================================")
  
}
# 
# 
current_folder = "/Volumes/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/"
inflow_list = list(RF04_inflow_summary, RF05_inflow_summary)
saveRDS(inflow_list, file.path(current_folder, "inflow_summary.RData"))
readRDS(file.path(current_folder, "inflow_summary.RData")) 
# a = list(RF04_n_pixel_summary, RF05_n_pixel_summary, RF05_all_conc_summary, RF04_all_conc_summary, RF05_small_conc_summary, RF04_small_conc_summary,  RF05_large_conc_summary, RF04_large_conc_summary)
# save(a, file = file.path(current_folder, "background_summary.RData"))
# #saveRDS(a, file = file.path(current_folder, "background_summary_RDS.RData"))
# #browser()
# saveRDS(RF04_n_pixel_summary, file.path(current_folder, "RF04_n_pixel_summary.RData"))
# saveRDS(RF05_n_pixel_summary, file.path(current_folder,"RF05_n_pixel_summary.RData"))
#sample_bg = raster::values(raster::crop(x1, e_current_background))
graphics.off()
# filename = paste("xch4", current.research.flight, "segment", current.segment.input, "nc",sep=".")
# filename.psurf = paste("psurface", current.research.flight, "segment",current.segment.input,"nc",sep=".")
# xy.facility= xy.release   ## nominal MiVida plant location
# sig.sel = en0   # sig mult above median for threshold
# nmin= n.sel    ## minimum number in clump
# xll=static.release.plot.lon
# yll=static.release.plot.lat          # larger plotting box
# inflow.box = static.release.inflow
# plume.box = static.release
# ccm.file = F
# l.plot = T
# plot.moniker = paste("RF05_release",  as.character(en0), as.character(segment.inputs[current.idx]))
#                                                
