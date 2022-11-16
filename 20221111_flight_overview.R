# IMPORT LIBRARIES
library(raster)
library(ggplot2)
library(reshape2)
library(ggthemes)
library(ggmap)
library(egg)
# library(scico)
library(gridExtra)
register_google(key = "AIzaSyBfk-8_ssaa77tiCB5G-a27QvthaIsAbcc")
#register_google(key = "AIzaSyBfk-8_ssaa77tiCB5G-a27QvthaIsAbcc", write = TRUE)
load("/Users/apisada/Documents/R_scripts/heat.275.RData")
cat('\014')
## LOCATION OF THE SOURCE (MiVida)
p_x = -103.465
p_y = 31.526
#########
rms.cal = function(input){
  output = sqrt(mean(input**2))
  return(output)
}
      
gplot_data <- function(x, maxpixels = 50000)  {
  x <- raster::sampleRegular(x, maxpixels, asRaster = TRUE)
  coords <- raster::xyFromCell(x, seq_len(raster::ncell(x)))
  ## Extract values
  dat <- utils::stack(as.data.frame(raster::getValues(x))) 
  names(dat) <- c('value', 'variable')
  
  dat <- dplyr::as.tbl(data.frame(coords, dat))
  
  if (!is.null(levels(x))) {
    dat <- dplyr::left_join(dat, levels(x)[[1]], 
                            by = c("value" = "ID"))
  }
  dat
}

###
# exported_map <- get_googlemap(c(p_x, p_y), maptype = "satellite", zoom = 5) 
# ggmap(exported_map)
# IMPORT DATA OF INTEREST 
current_folder = "/Users/apisada/Google Drive/My Drive/ON-GOING PAPER/Figures/"
raw_raster_folder = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/RF06/SEGMENT_RASTERS"
mosaic_raster_folder = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/RF06/Mosaic/"


# AIM FOR RF06 --> one of the MiVida overpass 
# CREAT PLOTS 
current_xch4 = raster(file.path(raw_raster_folder, "xch4.RF06.segment.3.nc")) *1e9
mosaic_xch4 = raster(file.path(mosaic_raster_folder, "RF06_mosaic.nc")) *1e9

mosaic_xch4_cropped = raster::crop(mosaic_xch4,  extent(c(-104.4, -102.8, 31.35, 32.55)))
plot(mosaic_xch4_cropped, col = heat.275)
# PLOT I
static.mivida = extent(c(-103.485, -103.455, 31.517, 31.545))
inflow.mivida = (raster::values(raster::crop(current_xch4, extent(c(-103.485, -103.455, 31.517, 31.520)))))

raw_enhanced_xch4 = raster::crop(current_xch4, static.mivida)
raw_enhanced_xch4_all = raw_enhanced_xch4

m0 = mean(inflow.mivida, na.rm =  T)
s0 = sd(inflow.mivida, na.rm =  T)
en0 = 1.5 

t0 = m0 + s0 * en0 

raster::values(raw_enhanced_xch4)[raster::values(raw_enhanced_xch4) < t0] = NA


raw_enhanced_clumped_xch4 = raw_enhanced_xch4 
raster::values(raw_enhanced_clumped_xch4) = NA

cl0 = clump( raw_enhanced_xch4 )
u0 = unique( raster::values(cl0) )
n.sel = 50 #pix
# loop over all clumps in the box, putting in xch4 for selected (>= 20 pix) clumps
for( u1 in u0[is.finite(u0)]) {
  l.sel = (raster::values( cl0 ) == u1 )
  l.sel[is.na(l.sel)] = F
  if( sum(raster::values( cl0 ) == u1, na.rm=T ) >= n.sel){
    plot(sum(raster::values( cl0 ) == u1, na.rm=T ))
    raster::values( raw_enhanced_clumped_xch4 )[ l.sel ] = raster::values( raw_enhanced_xch4 )[l.sel]
  }
}

Pts = c(p_x, p_y)
crop.radius = 1200 # 1 km
distanceFromPoints(raw_enhanced_clumped_xch4, Pts) < crop.radius -> cropped.area 
raw_enhanced_clumped_xch4[cropped.area == 0] <- NA

raw_enhanced_clumped_xch4 = raw_enhanced_clumped_xch4 - m0 


xch4_raw_df <- as.data.frame(raw_enhanced_xch4_all, xy = TRUE, na.rm = TRUE)

xch4_df <- as.data.frame(raw_enhanced_clumped_xch4, xy = TRUE, na.rm = TRUE)


mosaic_xch4_df <- as.data.frame(mosaic_xch4, xy = TRUE, na.rm = TRUE)

#mosaic_xch4_raw_gg <-ggplot( as.data.frame(mosaic_xch4, xy = TRUE, na.rm = TRUE), aes(x = x, y = y)) +

mosaic_xch4_raw_gg <-ggplot(mosaic_xch4_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = xch4_rect_filt_fill,)) + theme_minimal() + 
  ylab("Latitude") + xlab("Longitude") + labs(fill = "Overview XCH4 (ppb)")
# 
p_mosaic = mosaic_xch4_raw_gg   + theme(legend.position = "bottom") +  
  geom_point(data = source_df, shape = 24, fill = "plum1", size = 3) + 
  scale_fill_gradientn(colours = heat.275, limits = c(1850, 1920))
rm(mosaic_xch4_raw_gg)
p_mosaic

# p_mosaic = mosaic_xch4_raw_gg  +  scale_fill_gradientn(colours = heat.275)
# p_mosaic
# p_mosaic + theme(legend.position = "bottom") #+ theme_bw() 
# 
# xch4_raw_gg <-ggplot(xch4_raw_df, aes(x = x, y = y)) +
#   geom_raster(aes(fill = layer)) + theme_minimal() + 
#   ylab("Latitude") + xlab("Longitude") + labs(fill = "XCH4 (ppb)")
# # 
# # xch4_raw_gg <- ggmap(exported_map, base_layer = ggplot(xch4_raw_df, aes(x = x, y = y))) +
# #   geom_raster(aes(fill = layer)) + theme_minimal() + 
# #   ylab("Latitude") + xlab("Longitude") + labs(fill = "XCH4 (ppb)")
# 
# # REGULAR PLOT ATTEMPT 
# par(mfrow=c(1, 3))
# plot(crop(current_xch4, static.mivida), col = heat.275)
# plot(raw_enhanced_clumped_xch4, col = heat.275)
# 
# # GGPLOT ATTEMPT 
# 
# map_center_x = (extent(xch4_raw_df)[1] + extent(xch4_raw_df)[2])/2
# map_center_y = (extent(xch4_raw_df)[3] + extent(xch4_raw_df)[4])/2
# 
# exported_map <- get_googlemap(center = c(lon = map_center_x, lat = map_center_y), maptype = "satellite", zoom = 15) 
# 
# 
# exported_map2 <- get_googlemap(center = c(lon = map_center_x, lat = map_center_y), maptype = "satellite", zoom = 14) 
# 
source_df = data.frame(x = p_x, y = p_y)

xch4_gg <- ggplot(xch4_raw_df, aes(x = x, y = y)) +
  geom_raster(aes(fill = layer)) + theme_minimal() +
  geom_point(data = source_df, shape = 24, fill = "plum1", size = 3) + 
  ylab("Latitude") + xlab("Longitude") + labs(fill = "XCH4 (ppb)")
# 
# ####### RAW #####################
p1 = xch4_gg + theme_minimal() + theme(legend.position = "bottom") +  scale_fill_gradientn(colours = heat.275) 

ggarrange(p1 +  ylab("Latitude") + xlab("Longitude") + labs(fill = "XCH4 (ppb)"), 
          p_mosaic, nrow = 1, labels = c("A","B"))
# ########## IME ####################
# p2 = ggmap(exported_map, base_layer =  ggplot(xch4_df, aes(x = x, y = y))) + geom_raster(aes(fill = layer)) + coord_cartesian() +  scale_fill_gradientn(colours = heat.275) #+  theme_minimal() 
# p2
# 
# # ggmap(exported_map, base_layer =  ggplot(xch4_raw_df, aes(x = x, y = y))) +
# #   geom_raster(aes(fill = layer)) + coord_cartesian() +  scale_fill_gradientn(colours = heat.275) + 
# #   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.0007 * bsize, xmax= p_x + 0.00045 * bsize, ymin= p_y - 0.0004 * bsize, ymax= p_y + 0.0007 * bsize) + 
# #   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.0005 * bsize, xmax= p_x + 0.00035 * bsize, ymin= p_y - 0.0003 * bsize, ymax= p_y + 0.0005 * bsize) + 
# #   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.00025 * bsize, xmax= p_x + 0.00025 * bsize, ymin= p_y - 0.0002 * bsize, ymax= p_y + 0.00025 * bsize) + 
# #   theme_minimal() 
# 
# ######## DI PLOT ##################
# p_xmin = xmin(raw_enhanced_xch4_all)
# p_xmax = xmin(raw_enhanced_xch4_all)
# 
# p_ymin = xmin(raw_enhanced_xch4_all)
# p_ymax = xmin(raw_enhanced_xch4_all)
# 
# # p3 =  xch4_raw_gg  +  scale_fill_gradientn(colours = heat.275) + geom_rect(aes(xmin= p_x - 0.0005 * bsize, xmax= p_x + 0.0005 * bsize, ymin= p_y - 0.0005 * bsize, ymax= p_y + 0.0005 * bsize), alpha=0.02, fill="blue")
# # p3
# 
# # 
# # p3 =  xch4_raw_gg  +  
# #   scale_fill_gradientn(colours = heat.275) + geom_rect(aes(xmin= p_x - 0.0005 * bsize, xmax= p_x + 0.0005 * bsize, ymin= p_y - 0.0005 * bsize, ymax= p_y + 0.0005 * bsize), alpha=0.005, fill="grey") +
# #   geom_rect(aes(xmin= p_x - 0.00025 * bsize, xmax= p_x + 0.00025 * bsize, ymin= p_y - 0.00025 * bsize, ymax= p_y + 0.00025 * bsize), alpha=0.005, fill="red")
# # p3
# bsize = 15
# p3_old =  xch4_raw_gg  +
#   scale_fill_gradientn(colours = heat.275) + annotate("rect", color="white", alpha = 0, xmin= p_x - 0.0007 * bsize, xmax= p_x + 0.00045 * bsize, ymin= p_y - 0.0004 * bsize, ymax= p_y + 0.0007 * bsize) + 
#   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.0005 * bsize, xmax= p_x + 0.00035 * bsize, ymin= p_y - 0.0003 * bsize, ymax= p_y + 0.0005 * bsize) + 
#  annotate("rect", color="white", alpha = 0, xmin= p_x - 0.00025 * bsize, xmax= p_x + 0.00025 * bsize, ymin= p_y - 0.0002 * bsize, ymax= p_y + 0.00025 * bsize)
# 
# 
# #+geom_rect(aes(xmin= p_x - 0.00025 * bsize, xmax= p_x + 0.00025 * bsize, ymin= p_y - 0.00025 * bsize, ymax= p_y + 0.00025 * bsize), alpha=0.005, fill="red")
# p3 = ggmap(exported_map, base_layer =  ggplot(xch4_raw_df, aes(x = x, y = y))) +
#   geom_raster(aes(fill = layer)) + coord_cartesian() +  scale_fill_gradientn(colours = heat.275) + 
#   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.0105, xmax= p_x + 0.00475, ymin= p_y - 0.004, ymax= p_y + 0.0105) + 
#   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.0075, xmax= p_x + 0.00425, ymin= p_y - 0.0035, ymax= p_y + 0.0075) + 
#   annotate("rect", color="white", alpha = 0, xmin= p_x - 0.00375, xmax= p_x + 0.00375, ymin= p_y - 0.003, ymax= p_y + 0.00375) + 
#   theme_minimal() 
# 
# p3
# 
# ##################################
# #-103.465   31.526
# grid.arrange(p1, p2 +  ylab("Latitude") + xlab("Longitude") + labs(fill = "XCH4 (ppb)"), p3 +  ylab("Latitude") + xlab("Longitude") + labs(fill = "XCH4 (ppb)"), nrow = 1)
# 
# ##################################
# # exported_map <- get_googlemap(center = c(lon = p_x, lat = p_y), maptype = "satellite", zoom = 14) 
# # #ggmap(exported_map)
# # ggmap(exported_map)
# # 
# # a = ggmap(exported_map) + annotate("rect", color="red", fill = "blue", alpha = 0.5, xmin= p_x - 0.00375, xmax= p_x + 0.00375, ymin= p_y - 0.00375, ymax= p_y + 0.00375) 
# # a
# # a = ggmap(exported_map, base_layer =  ggplot(xch4_raw_df, aes(x = x, y = y))) +
# #   geom_raster(aes(fill = layer)) + coord_cartesian() +  scale_fill_gradientn(colours = heat.275) + annotate("rect", color="white", fill = "blue", alpha = 0, xmin= p_x - 0.00375, xmax= p_x + 0.00375, ymin= p_y - 0.00375, ymax= p_y + 0.00375) + theme_minimal() 
# # a
# 
# 
# p_mosaic = ggmap(exported_map2, base_layer =  ggplot(mosaic_xch4_df, aes(x = x, y = y))) + geom_raster(aes(fill = layer)) + coord_cartesian() +  scale_fill_gradientn(colours = heat.275) #+  theme_minimal() 
# p_mosaic
# 
# # Have all 