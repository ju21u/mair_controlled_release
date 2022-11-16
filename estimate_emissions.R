rms.cal = function(input){
  output = sqrt(mean(input**2))
  return(output)
}

estimate_emissions = function( current.research.flight = "RF04", xy.input, box.inflow.input, 
                               e.input, source.name.input, current.segment.input, res = "5x1"){
  
current.segments = current.segment.input
Pts = matrix( unlist(xy.input), ncol=2, nrow=1, byrow =T)
e.inflow = extent (box.inflow.input)
e.current = e.input
source.name = source.name.input
#source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF04/load_functions_constants.R")

if (current.research.flight == "RF04"){
  current.date = "2021-07-30"
  if (res == "1x1"){
  current.path = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF04/1x1/SEGMENT_RASTERS/"
  } else if (res == "5x1"){
    current.path = "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF04/5x1/SEGMENT_RASTERS/"
  }
  source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/Bella/20211028_load_data_RF04.R")

}else if (current.research.flight == "RF05"){
  current.date = "2021-08-03"
   source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF05/20211101_load_data_temp.R")
  if (res == "1x1"){
    current.path =  "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF05/1x1/SEGMENT_RASTERS"
  } else if (res == "5x1"){
    current.path =  "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/all.segs/RF05/5x1/SEGMENT_RASTERS"
  }
  # source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF04/calculate_les.R")
}else if (current.research.flight == "RF06"){
  current.date = "2021-08-06"
  source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF06/20220213_load_data_RF06.R")
  if (res == "1x1"){
    current.path =   "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/MiVida/RF06/"
  } 
  
}else if (current.research.flight == "RF07"){
  current.date = "2021-08-09"
  source("/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Scripts/clean_RF07/20220213_load_data_RF07.R")
  if (res == "1x1"){
    current.path =  "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/MiVida/RF07/"
  } 
  
}



  
  
  #################
  #################
  
  ## find the medians, std dev, of the inflow in the release box, for each segment.  
  ## then define the threshold ("thr") for image identification/deåinition  for each 
  ##par(mfrow =c( 1,2))å
  my.median = NULL ; my.sd = NULL ; my.thr = NULL # initialize
  my.clumps = NULL

  
  for  ( k in current.segments ){ #, 11, 13, 15, 17, 19, 21, 23, 25, 3, 30, 5, 7, 9 )  ) { #c(1, 5, 10, 12, 14) ) {
    # open the png, the results will be committed to this figure
    # read th xch4, psurface (needed for mass computation), and albedo rasters for the segment
   # browser()
    setwd(current.path)
    #plot(stack("regridded_diag_0p85pixwidth_vida-pass0.nc",varname="wt"))
    #r = raster::raster("regridded_diag_0p85pixwidth_vida-pass0.nc")
    #raster::crs(r) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
    x1 = raster(paste("xch4", current.research.flight, "segment", k, "nc",sep="."))
    ## crop them to the release box
    x1e = raster::crop( x1, e.current)
    p1e = raster::crop( raster(paste("psurface", current.research.flight, "segment",k,"nc",sep=".")), e.current)
    #calculate the mean of the entire domain?
    #browser()
    p1.mean = rms.cal(raster::values(p1e)[!is.na(raster::values(p1e))]) #rms.cal( raster::values(p1e)) # , na.rm=T)
    #p1.mean = 920
    assign(paste("release.seg", k, current.research.flight, sep="."), x1e )
    ## cropt to the inflow box
    x2 = raster::crop(x1, e.inflow)
    m0 = median( raster::values(x2), na.rm=T ) 
    s0 =  sd( raster::values(x2), na.rm=T ) 
   # browser()
    ## ----  median and sd for the inflow box
    # t0 is the threshold for maksing --> plume identification and IME
    t0 = m0 + s0 * en0 # maskingraster::values
    my.median = c(my.median, m0) 
    my.sd = c(my.sd, s0) 
    my.thr = c(my.thr,  t0)
    ## x1m is the masked box -- shows only the plume data, other pixels = NA
    x1m = x1e
    raster::values(x1m)[raster::values(x1e) < t0] = NA
    # cl0 is the raster of clumps, each pixel with value identifiying the clump number (ID) or NA if not in a clump
    
    
    x1m.3 <- x1m
    # CROP RADIUS
    # CROP RADIUS
    distanceFromPoints(x1m.3, Pts) < crop.radius -> cropped.area # SELECTED
    # CROP RADIUS
    # CROP RADIUS
    #()
    x1m.3[cropped.area == 0] <- NA
    #plot(x1m, main = "original")
    plot(x1m.3, main = paste(source.name.input, "cropped"), col = heat.275)
    plot(x1m, main = paste("Segment", k, "(original,", current.research.flight, source.name.input, ")"), col = heat.275)
    plot(x1m.3, main =  paste("Segment", k, "(cropped", current.research.flight, source.name.input, ")"), col = heat.275 )		  
    x1m <- x1m.3
    plot(x2, main = 'obs inflow')
    x1m.2 = x1m  ; raster::values( x1m.2) = NA
    
    cl0 = clump( x1m )
    u0 = unique( raster::values(cl0) )
    cl2 = cl0    ; raster::values( cl2 ) = NA
    # loop over all clumps in the box, putting in xch4 for selected (>= 20 pix) clumps
    for( u1 in u0[is.finite(u0)]) {
      l.sel= (raster::values( cl0 ) == u1 )
      l.sel[is.na(l.sel)] = F
      if( sum(raster::values( cl0 ) == u1, na.rm=T ) >= n.sel) raster::values( x1m.2 )[ l.sel ] = raster::values( x1m )[l.sel]
      if( sum(raster::values( cl0 ) == u1, na.rm=T ) >= n.sel) raster::values( cl2 )[ l.sel ] = raster::values( cl0 )[l.sel]
    }
    x1m = x1m.2  ## x1m is the masked xch4 clumps, now selected for the number of points >= n.sel 
    #browser()
    
    png(file = file.path("/Volumes/GoogleDrive/My Drive/ON-GOING PAPER/Figures/plumes/", paste("cropped", current.research.flight, source.name.input, paste(k, ".png", sep =""), sep = "_") ))
    plot(x1m.2 * 1e9, main =  paste( current.research.flight, "Segment", k, source.name.input), col = heat.275,
         legend.args=list(text='XCH4 (ppb)', side=4, font=2, line=2.8, cex=0.8))		
    
#    points(xy.release$x, xy.release$y, pch = 21, col = "red")
    points(xy.release$x, xy.release$y, pch = 17, col = "grey", cex = 1.2)
    dev.off()
    #browser()
    cl0 = cl2   ## cl0 is the masked clump ID wraster, now selected for the number of points >= n.sel 
    u0 = unique(raster::values( cl0))
    #browser() 
    #n0 = NULL
    # select each plume in the box, comput IME, area, etc
    for( u1 in u0[is.finite(u0)]) {
      #n0 = rbind( n0, c(u1, sum(rvalues( cl0 ) == u1, na.rm=T ) ) ) 
      l.sel= (raster::values( cl0 ) == u1 )
      l.sel[is.na(l.sel)] = F
      xy0 = coordinates (cl0 )
      # center of the plume area
      mean.lat= rms.cal(xy0[l.sel,"y"]) #, na.rm=T)
      mean.lon= rms.cal(xy0[l.sel,"x"]) #, na.rm=T)
      # sum(raster::values(raster::area(x1m))) * 1e6 # 
      a0 = sum(raster::values(raster::area(cl0))[l.sel], na.rm=T) * 1e6
      # distance of the plume *center* to the release point
      dist0 = spDistsN1( pt=c(mean.lon, mean.lat), pts =Pts, longlat=T) * 1e3
      kg = sum( molef.to.kg ( (raster::values(x1m)[l.sel] - m0 ), p1.mean) *raster::values(raster::area( x1m))[l.sel] * 1e6, na.rm=T)
      signals = mean(raster::values(x1m)[l.sel] - m0 )
      noise = cellStats(x2, sd)
      l.sel.n = sum(l.sel) #()
      res = c("seg."=k, "clump"=u1, "lat"=mean.lat,"lon"=mean.lon,  "IME.kg"=kg, 
              "dist.m"=dist0, "Tot.area"=a0, "infL.med"=m0, "infL.sd"=s0, 
              "Sel.thr"=t0, "signal" = signals, "noise" = noise, "l.sel.n" = l.sel.n)
      if(kg > 0){
        my.clumps = rbind( my.clumps, res ) }
    }
    print(paste('FINISHED 1 SEGMENT', as.character(k)))
  }
  #write.table(signif(my.clumps,6), file="/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/Intermediates/summary_csvs/RF07_summary.csv",row.names=F,col.names=T,quote=F)
  #setwd( old.dir )
  
  #browser()
  # ======================================================================================
  #cat('\014')
  nc_path =  "/Users/apisada/Google Drive/My Drive/Research/Harvard_Research/MethaneAIR/"
  if (!is.null(my.clumps)){
    summary.current.source = signif(my.clumps,6)
    View(summary.current.source)
    #xy.basking = list("x" =  -103.696 , "y" =  32.366)
    
    # box.basking  = c(-103.54, -103.47, 31.557, 31.515)
    # e.basking.prime = extent (box.basking)
    # 
    # box.inflow.basking.prime = c(-103.54, -103.47, 31.557, 31.566)
    # 
    # e.inflow = extent (box.inflow.basking.prime)
    # e.current = e.basking.prime
    
    
    summary.current.RF <- read.csv(file.path(nc_path, "From_Steve/all.segs/", current.research.flight, "/5x1/all.segs.csv"), sep = " ")
    summary.current.RF <- summary.current.RF[rowSums(is.na(summary.current.RF)) == 0, ]   
    if (current.research.flight == "RF05"){
      summary.current.RF$seg <- seq(length(summary.current.RF$start))}
    #========================================
    df_summary <- data.frame(matrix(vector(), 0, 13,
                                    dimnames=list(c(), c("SEGMENT", "IME EMISS", "IME AREA", "IME IME", "clumps", "IME300 EMISS", "IME300 AREA", "IME300 IME",  "clumps",  "IME COEF EMISS", "IME 300 COEF EMISS", " WIND", "THRES"))),
                             stringsAsFactors=F)
    df_summary_condense <- data.frame(matrix(vector(), 0, 17,
                                             dimnames=list(c(), c("SEGMENT", "IME EMISS", "IME AREA", "IME IME", "clumps", 
                                                                  "RATIO EMISS", "RATIO SD", "RATIO 25",  "RATIO 975", 
                                                                  "IME 25", "IME 975", "SD"," IME COEF EMISS", 
                                                                  "IME 300 COEF EMISS"," WIND", "THRES", "SNR"))),
                                      stringsAsFactors=F)
    #========================================
    seg.numbers = current.segments# c(1) #, 5, 10, 12)
    RF.df <- data.frame(seg.numbers= seg.numbers,
                        date = rep(current.date, each = length(seg.numbers)),
                        t.start =as.numeric(summary.current.RF$start[is.element(summary.current.RF$seg, seg.numbers)]),
                        t.stop =as.numeric(summary.current.RF$stop[is.element(summary.current.RF$seg, seg.numbers)]))
    
    current.thres = en0
    
    
    ###########   
    #browser()
    keep.idx = c()
    for (seg.idx in seq(length(seg.numbers))){
      print('starting a new loop')
      seg.number = RF.df$seg.numbers[seg.idx]
      t.start =  as.POSIXct(strptime(paste(RF.df$date[seg.idx], floor(RF.df$t.start[seg.idx]), floor((RF.df$t.start[seg.idx] %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")
      t.start.idx = which.min(abs(Times - t.start)) - 15
      t.stop =  as.POSIXct(strptime(paste(RF.df$date[seg.idx], floor(RF.df$t.stop[seg.idx]), floor((RF.df$t.stop[seg.idx] %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")
      #t.stop =  as.POSIXct(strptime("2021-07-30_18:00:00", "%Y-%m-%d_%H:%M:%S"), tz = "UTC")
      t.stop.idx = which.min(abs(Times - t.stop)) + 15 
      if (t.stop.idx > length(Times)){
        t.stop.idx = length(Times)
      }
      # if (t.start.idx == t.stop.idx){
      #   t.stop.idx = t.stop.idx + 1
      # }
      print('done selecting the time')
      #browser()
      rr.stack.crop <- rr.stack #raster::crop(rr.stack, e.current)
      rr.t <- subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5))
      #browser()
      rr.t.mean <- raster::calc( subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = rms.cal)
      
      #rr.t.mean <- calc( subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = mean, na.rm = T)
      # if (dim(subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)))[3] > 1){
      #   rr.t.mean <- raster::calc( subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = mean, na.rm = T)
      # } else if (dim(subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)))[3] == 1){
      #   rr.t.mean <- rr.stack.crop
      # } else{print("there is sth wrong")}
      
      plot(rr.t.mean, col = heat.275, main = paste("LES: Segment", seg.number, source.name, current.research.flight))
      
      
      rr.inflow.stack.crop <- raster::crop(rr.stack, e.inflow)
      rr.inflow.t.mean <- raster::calc( subset(rr.inflow.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = rms.cal)
      
      psf.inflow.stack.crop <- raster::crop(psfc.stack, e.inflow)
      psd.inflow.t.mean <-  subset(psf.inflow.stack.crop,  seq(t.start.idx, t.stop.idx, 5))
      psd.inflow.t <- raster::calc( subset(psf.inflow.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = rms.cal)
      
      n_matched = sum(summary.current.source[, 'Tot.area'], na.rm = T)/median(raster::values(raster::area(rr.t.mean)) * 1e6)
      #browser()
      l.sel= (raster::values( rr.t.mean ) > sort(raster::values(rr.t.mean) )[1])
      l.sel[is.na(l.sel)] = F
      les.mass = sum( molef.to.kg ( (sort(raster::values(rr.t.mean) )[1:n_matched] ),  mean(raster::values(psd.inflow.t.mean), na.rm =  T)) ) * 0.111 * 0.111 # * raster::values(raster::area(rr.t.mean))[l.sel] * 1e6, na.rm=T)
      # sum( molef.to.kg ( (raster::values(x1m)[l.sel] - m0 ), p1.mean) *raster::values(raster::area( x1m))[l.sel] * 1e6, na.rm=T)
      les.mass.all = NULL
      #browser()
      #graphics.off()
      for (idx_layer in seq(1, nlayers(rr.t))){
        #print(n_matched)
        x = subset(rr.t, names(rr.t)[idx_layer] )
        raster::plot(x)
        #contour(x, maxpixels = n_matched,  add = TRUE)
        if (exists("x") ){
          if(n_matched < 1 & n_matched > 0){
            n_matched = 1
          }
            
          selected_contourline <- rasterToContour(x, levels = sort(raster::values(x) , decreasing = T)[n_matched] )
          plot(selected_contourline, col = '#f78605', lwd = 2, add = TRUE)
          }
        les.mass.all = c(les.mass.all, sum( molef.to.kg ( (sort(raster::values(x) , decreasing = T)[1:n_matched] ),  mean(raster::values(psd.inflow.t), na.rm =  T))* 0.111 * 0.111, na.rm=T))
      }
     # browser()
      plot(sort(raster::values(x1m), decreasing = TRUE), sort(raster::values(x), decreasing = TRUE) [1: length(sort(raster::values(x1m), decreasing = TRUE))],
           xlab = "obs.xch4", ylab = "les.ch4", main = "obs vs les")
      
        #molef.to.kg(sum(sort(raster::values(rr.t.mean), decreasing = T)[1:n_matched], na.rm =  T), mean(raster::values(psd.inflow.t.mean), na.rm =  T))
     #kg = sum( molef.to.kg ( (raster::values(x1m)[l.sel] - m0 ), p1.mean) *raster::values(raster::area( x1m))[l.sel] * 1e6, na.rm=T)
      obs.mass =  sum(summary.current.source[, 'IME.kg'], na.rm = T)
      #ratio.estimate = obs.mass/les.mass * 988
      ratio.estimate.all = obs.mass/les.mass.all * 988
      #ratio.estimate = rms.cal(ratio.estimate.all)
      ratio.sd =  sd(ratio.estimate.all)
      #qs.emisss = quantile(ratio.estimate.all, c(0.25, 0.75), na.rm = TRUE)
      bt.mean = bootstrap::bootstrap(ratio.estimate.all, 1000, rms.cal)$thetastar
      
      
      
      #(bt.mean, 0.025) #
      #browser()
      ratio.estimate = quantile(bt.mean, 0.5) 
      ratio25 = quantile(bt.mean, 0.025) #qs.emisss[1]
      ratio975 = quantile(bt.mean, 0.975) #qs.emisss[2]
      
      bt.prc = bootstrap::bootstrap(sqrt(((ratio.estimate.all - ratio.estimate)/ratio.estimate)**2)*100, 1000, rms.cal)$thetastar
      #browser()
      # print("========== Look =======")
      # print(ratio.estimate)
      #browser()
      #+++++++++++++++++++++++++++++++
  
      
      #+++++++++++++++++++++++++++++++
      # if (dim(subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)))[3] > 1){
      # rr.t.mean <- raster::calc( subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = mean, na.rm = T)
      # } else if (dim(subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)))[3] == 1){
      #   rr.t.mean <- rr.stack.crop
      # } else{print("there is sth wrong")}
      
      plot(rr.inflow.t.mean, col = heat.275, main = paste("Segment", seg.number, "(inflow", source.name.input, ")"))
     # browser()
      
      #u.stack.crop <- raster::crop(u.stack, e.inflow)
      #u.t.mean <- calc( subset(u.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = mean, na.rm = T)
      #plot(u.t.mean)
      
      #v.stack.crop <- raster::crop(v.stack, e.inflow)
      wind.stack.crop <- raster::crop(wind.stack, e.inflow)
      wind.stack.crop.t <-  subset(wind.stack.crop,  seq(t.start.idx, t.stop.idx, 5))
      #wind.t.mean <- calc( wind.stack.crop.t, fun = rms.cal)
      
      
      
      # IME 
      current.wind =  rms.cal(raster::values(wind.stack.crop.t)) #cellStats(wind.t.mean, "rms")
      current.winds = raster::values(wind.stack.crop.t) #cellStats(wind.stack.crop.t, "rms")
      #current.winds = cellStats(wind.stack.crop.t, "rms") 
      
      current.ime =  sum(summary.current.source[, 5][summary.current.source[, 1] ==  seg.number])
      current.imes =  summary.current.source[, 5][summary.current.source[, 1] ==  seg.number]
      current.areas  = summary.current.source[, 7][summary.current.source[, 1] ==  seg.number & summary.current.source[, 5] > 0]
      current.area  = sum(summary.current.source[, 7][summary.current.source[, 1] ==  seg.number & summary.current.source[, 5] > 0])
      print("Calculating the emissions and quantiles")
      #ime.emiss = emission.cal.IME(wind = current.wind, ime_input = current.ime, area_input = current.area)
      ime.emisss = emission.cal.IME(wind = current.winds, ime_input = current.ime, area_input = current.area)
      #ime.emiss = rms.cal(ime.emisss)
      bt.mean = bootstrap::bootstrap(ime.emisss, 1000, rms.cal)$thetastar
      #(bt.mean, 0.025) #
      ime.emiss = quantile(bt.mean, 0.5)
      min25 = quantile(bt.mean, 0.025) #qs.emisss[1]
      max975 = quantile(bt.mean, 0.975) #qs.emisss[2]
      sd5 =  sd(ime.emisss)
      print(ime.emiss)
      #browser()
      # coef.ime.emiss = emission.cal.IME(wind = current.wind, ime_input = current.ime, area_input = current.area,  coef = 0.3349, interp = 0.7521)
      # print(coef.ime.emiss)
      
      #coef.ime.emiss = emission.cal.coef.IME(wind = current.wind, ime_input = current.ime,  area_input = current.area, les_ime_input = les.ime.input, les_area_input = les.area.input, les_emiss_hr = 988)
      #print("new coef 300 (wrong):")
      #print(coef.ime.emiss)
      
      print("===================================")
      mean.signal = mean(summary.current.source[, 'signal'])
      #weighted mean
      sd.noise = mean(summary.current.source[, 'noise'] )
      SNR = mean.signal/sd.noise
      # Ratio from Steve 
      #	names(rr.stack) <- seq(421)
      # plot(crop(subset(rr.stack, "X200"), e))
      a.list = list(seg.number, ime.emiss,  current.areas, current.imes, length(current.imes),  ime.emiss, current.areas, current.imes, length(current.imes), 0, 0, current.wind, current.thres, SNR)
      if (sum(is.na(a.list)) == 0){
        df_row = data.frame(a.list)
        names(df_row)<- names(df_summary)
        df_summary <- rbind(df_summary, df_row)
        
        df_row = data.frame(list(seg.number, ime.emiss,  current.area, current.ime, length(current.imes),  ratio.estimate, ratio.sd, ratio25, ratio975, min25, max975, sd5, 0, 0, current.wind, current.thres, SNR))
        names(df_row)<- names(df_summary_condense)
        df_summary_condense <- rbind(df_summary_condense, df_row)
        
        
        rt <- rts(rr.stack, Times) # creating a RasterStackTS object
        #plot(subset(rt, seq(t.start.idx, t.stop.idx, 5)), col = heat.275)
        keep.idx = c(keep.idx, seg.idx)
      }else (print('skip'))
      # df_summary <- data.frame(matrix(vector(), 0, 9,
      #                                 dimnames=list(c(), c("SEGMENT", "IME EMISS", "IME AREA", "IME IME", "IME300 EMISS", "IME300 AREA", "IME300 IME", " WIND", "THRES"))),
      #                          stringsAsFactors=F)
      print('finished one loop')
    }
    print(paste("Done with this execution: ", current.research.flight, source.name, as.character(current.segment.input)))
    print("==============================================================")
    df_summary_condense %>% mutate_at(vars(IME.EMISS, IME.AREA, IME.IME, RATIO.EMISS, RATIO.SD, RATIO.25, RATIO.975, IME.25, IME.975, SD, X.WIND), function(x)(round(x, 2))) -> df_summary_condense
    
    rownames(df_summary_condense) <- seq(dim(df_summary_condense)[1])
    df_summary_condense$formatted.time  <- as.POSIXct(strptime(paste(RF.df$date, floor(RF.df$t.start), floor((RF.df$t.start %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")[keep.idx]
    #write.csv(df_summary_condense,file.path("/Volum∂ç∂es/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/Outputs/CSVs/", paste('df_summary_', current.research.flight, '_', source.name, '_thres_', toString(current.thres), '.csv', sep = '' )))
    RF.df$formatted.time <- as.POSIXct(strptime(paste(RF.df$date, floor(RF.df$t.start), floor((RF.df$t.start %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")
    View(df_summary_condense)
}
  else{df_summary_condense = NULL }
  
  #browser()
  return(list(df_summary_condense,  cbind(les.mass.all, as.vector(current.winds)) ))
}

