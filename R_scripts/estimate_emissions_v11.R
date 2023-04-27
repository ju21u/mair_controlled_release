estimate_emissions = function( current.research.flight = "RF04", 
                               xy.input, box.inflow.input, 
                               e.input, source.name.input, 
                               current.segment.input, 
                               res = "5x1", N_i = 100, 
                               plume_crop = NULL, raster_input = NULL, 
                               pressure_input = NULL, skip_winds = NULL,
                               resample_estimates = NULL,
                               unit_ppb = FALSE, DI_box = FALSE, 
                               save_raw = FALSE){
    
    # STEP 1 #######
    ## Convert the input of the functions to usable format within the function
    Pts = matrix( unlist(xy.input), ncol=2, nrow=1, byrow =T) 
    # extent of the inflow (excluding the plume)
    e.inflow = extent (box.inflow.input)
    e.current = e.input # extent of the current plume (excluding the inflow)
    source.name = source.name.input
    # if(sum(raw_maryann_data$index == current.idx) > 0){
    #   e.original = e.current
    #   #current_box = unlist(raw_maryann_data[raw_maryann_data$index == current.idx, 9:12])
    #   e.current = extent(unlist(raw_maryann_data[raw_maryann_data$index == current.idx, 9:12]))
    #   print("Use Maryann's extent!")
    #   #unlist(raw_maryann_data[raw_maryann_data$index == current.idx, 9:12])
    # }
   
    #################
    #################
    
    
    # STEP 2 #######
    # PREPARING DATA FOR RESAMPLING
    ## find the medians, std dev, of the inflow in the release box, 
    # for each segment.  
    ## then define the threshold ("thr") for image identification
   # browser()
    setwd(current.path)
    if (is.null(raster_input)){
    entire_raster_input = raster(paste("xch4.", current.research.flight,  
                                       ".segment.", current.segment.input, 
                                       ".nc",sep=""))
    } else{entire_raster_input = raster(raster_input)}
   
   # browser()
    if(unit_ppb){
      entire_raster_input = entire_raster_input/(1e9)
    }

    source_raster = raster::crop( entire_raster_input, e.current)
    
    
    if (is.null(pressure_input)){
    entire_input_pressure =  raster(paste("psurface.", current.research.flight,  ".segment.", current.segment.input, ".nc",sep=""))
    }else{entire_input_pressure = raster(pressure_input) }
    
   
    cropped_input_pressure = raster::crop(entire_input_pressure, e.current)
    entire_pressure_input = rms.cal(raster::values(cropped_input_pressure)[!is.na(raster::values(cropped_input_pressure))]) + rnorm(1, mean = 0, sd = sd(raster::values(cropped_input_pressure)[!is.na(raster::values(cropped_input_pressure))], na.rm = T))
    
    assign(paste("release.seg", current.segment.input, current.research.flight, sep="."), source_raster )
    #plume_shape = raster::raster(file.path(cropped_plume_folder, paste("cropped_", current.research.flight, "_release_", current.segment.input, ".nc", sep = "" )))
   
    
    ## crop to the inflow box
    
    source_inflow = raster::crop(entire_raster_input, e.inflow) 
    plot(source_inflow, main = "inflow", col = heat.275)
    # NEED TO MCMC THIS LINE ABOVE ^^^ !!!!!
    inflow_median = median( raster::values(source_inflow), na.rm=T ) 
    inflow_sd =  sd( raster::values(source_inflow), na.rm=T ) 
    
    SOURCE_INFLOW_MEDIAN = median(raster::values(source_inflow), na.rm=T) 
    SOURCE_INFLOW_SD =  sd(raster::values(source_inflow), na.rm=T) 
    # current sd of the original selected scene
    
    original_threshold = SOURCE_INFLOW_MEDIAN + SOURCE_INFLOW_SD * current_thres
    
    counter = 0 # counter of the iterations
    
    source_raster_save = source_raster 
#    plot(crop(entire_raster_input, bigbox), col = heat.275,
#         main = paste(current.research.flight, current.segment.input,
#                      "mIME extent only"))
    rect(e.current[1],  e.current[3], e.current[2],
         e.current[4], border = "blue")
    points(Pts, col = "red", pch = 15)
    
    # make a copy of the original source raster
#    plot(crop(entire_raster_input, bigbox), col = heat.275,
#         main = paste(current.research.flight, current.segment.input,
#                      "DI & mIME extents"))
    if(DI_box){
    current_box = unlist(raw_maryann_data[raw_maryann_data$index == current.idx, 9:12])
    rect(current_box[1],  current_box[3], current_box[2],
         current_box[4], border = "green")
    }
    rect(e.current[1],  e.current[3], e.current[2],
         e.current[4], border = "blue", lty = 2)
    
        # plot(source_raster, col = heat.275, main = "raw")
    points(Pts, col = "red", pch = 15)
    #browser()
    # STEP 3 #######
    # INITIALIZE THE OUTPUTS
    # OUTPUT DATAFRAME
    df_summary_combined <- data.frame(matrix(vector(), 0, 18,
                          dimnames=list(c(), c("SEGMENT", "IME EMISS", 
                          "IME AREA", "IME IME", "clumps", 
                          "RATIO EMISS", "RATIO SD", "RATIO 25",  "RATIO 975", 
                          "IME 25", "IME 975", "SD"," IME COEF EMISS", 
                          "IME 300 COEF EMISS"," WIND", "THRES", "SNR", 
                          "formatted.time"))),
                                    stringsAsFactors=F)
    # OUTPUT VECTORS
    my.median = NULL ; my.sd = NULL ; my.thr = NULL ; selected_clumps = NULL
    les.mass.all = NULL
    current.winds = NULL
    raw_row = NULL
    # STEP 4: APPLYING THRESHOLD, CROPPING THE SOURCE, CLUMPING ALGORITHM, and FILTERING SMALL CLUMPS 
    source_above_thres = source_raster
    # >> THRESHOLDING 
    # >> THRESHOLDING 
    raster::values(source_above_thres)[raster::values(source_raster) < original_threshold] = NA
    # >> THRESHOLDING 
    # >> THRESHOLDING 
 #   plot(source_raster, col = heat.275, main = "thres only")
    temp_source_above_thres_radius = source_above_thres
    # CROP AROUND THE SOURCE 
    # CROP AROUND THE SOURCE 
    
    # CROP AROUND THE SOUR CE + 0.001 degree 

    # CROP AROUND THE SOUR CE + 0.001 degree 
    
    # CROP EXACTLY AROUND THE SOURCE
    #temp_resample_source_above_thres_radius <- resample_source_above_thres * (plume_shape > 0)
    # CROP EXACTLY AROUND THE SOURCE
    #plot(temp_resample_source_above_thres_radius)
    
    # CROP AROUND THE SOURCE 
    # CROP AROUND THE SOURCE 
    
    # >> CROP RADIUS
    # >> CROP RADIUS
    distanceFromPoints(temp_source_above_thres_radius, Pts) < crop.radius -> cropped.area 
    
    # >> CROP RADIUS
    # >> CROP RADIUS
    
    
    
    temp_source_above_thres_radius[cropped.area == 0] <- NA
    # NA those pixels outside the radius 
    source_above_thres = temp_source_above_thres_radius
    
    temp_source_above_thres_clump = temp_source_above_thres_radius  
    raster::values(temp_source_above_thres_clump) = NA
    # >> CLUMPING
    # >> CLUMPING
    all_clump_IDs = clump(source_above_thres)
    # >> CLUMPING
    # >> CLUMPING
    unique_clump_IDs = unique( raster::values(all_clump_IDs) )
    temp_all_clump_IDs = all_clump_IDs    
    raster::values( temp_all_clump_IDs ) = NA

    # loop over all clumps in the box, putting in xch4 for selected (>= 20 pix) clumps
    for( current_unique_clump_ID in unique_clump_IDs[is.finite(unique_clump_IDs)]) {
      
      selected_clump_pixels= (raster::values( all_clump_IDs ) == current_unique_clump_ID )
      selected_clump_pixels[is.na(selected_clump_pixels)] = F
      
      if( sum(raster::values( all_clump_IDs ) == current_unique_clump_ID, na.rm=T ) >= min_pixels_in_a_clump) {
        raster::values( temp_source_above_thres_clump )[ selected_clump_pixels ] = raster::values( source_above_thres )[selected_clump_pixels]
        raster::values( temp_all_clump_IDs )[ selected_clump_pixels ] = raster::values( all_clump_IDs )[selected_clump_pixels]}
    }
    source_above_thres = temp_source_above_thres_clump  ## resample_source_above_thres is the masked xch4 clumps, now selected for the number of points >= min_pixels_in_a_clump 
    all_clump_IDs = temp_all_clump_IDs   ## all_clump_IDs is the masked clump ID wraster, now selected for the number of points >= min_pixels_in_a_clump 
    unique_clump_IDs = unique(raster::values( all_clump_IDs))
    
    
    plot(source_above_thres * 1e9, col = heat.275, 
         main = paste(current.research.flight, " Plume", current.idx, "(ppb)"))
    points(Pts, col = "red", pch = 15)
   browser()
    for( current_unique_clump_ID in unique_clump_IDs[is.finite(unique_clump_IDs)]) {
      selected_clump_pixels = (raster::values( all_clump_IDs ) == current_unique_clump_ID )
      selected_clump_pixels[is.na(selected_clump_pixels)] = F
      xy0 = coordinates (all_clump_IDs )
      
      # center of the plume area
      mean.lat= rms.cal(xy0[selected_clump_pixels,"y"]) #, na.rm=T)
      mean.lon= rms.cal(xy0[selected_clump_pixels,"x"]) #, na.rm=T)
      
      area_of_the_clump = sum(raster::values(raster::area(all_clump_IDs))[selected_clump_pixels], na.rm=T) * 1e6
      # distance of the plume *center* to the release point
      
      dist0 = spDistsN1( pt=c(mean.lon, mean.lat), pts =Pts, longlat=T) * 1e3
      kg = sum( molef.to.kg ( (raster::values(source_above_thres)[selected_clump_pixels] - inflow_median ), entire_pressure_input) *raster::values(raster::area( source_above_thres))[selected_clump_pixels] * 1e6, na.rm=T)
      signals = mean(raster::values(source_above_thres)[selected_clump_pixels] - inflow_median )
      noise = cellStats(source_inflow, sd)
      number_of_selected_clump_pixels = sum(selected_clump_pixels) #()
      
      current_selected_clump = c("seg."= current.segment.input, "clump"=current_unique_clump_ID, "lat"=mean.lat,"lon"=mean.lon,  "IME.kg"=kg, 
                                 "dist.m"=dist0, "Tot.area"=area_of_the_clump, "infL.med"=inflow_median, "infL.sd"=inflow_sd, 
                                 "Sel.thr"=original_threshold, "signal" = signals, "noise" = noise, "selected_clump.n" = number_of_selected_clump_pixels)
      
      
      if(kg > 0){
        selected_clumps = rbind( selected_clumps, current_selected_clump ) }
    }
    skipping = TRUE

    
   
    #
    # STEP 5 ##### 
    ## COMBINING ALL THE SELECTED CLUMPS 

    if (!is.null(selected_clumps)){
      skipping = FALSE
      current_row_df_summary <- data.frame(matrix(vector(), 0, 17,
                                                  dimnames=list(c(), c("SEGMENT", "IME EMISS", 
                                                                       "IME AREA", "IME IME", "clumps", 
                                                                       "RATIO EMISS", "RATIO SD", "RATIO 25",  
                                                                       "RATIO 975", "IME 25", "IME 975", "SD",
                                                                       " IME COEF EMISS", "IME 300 COEF EMISS",
                                                                       " WIND", "THRES", "SNR"))),
                                           stringsAsFactors=F)
      
      
      summary.current.source = signif(selected_clumps,6)
      
      #========================================
    # browser()
    if(is.null(skip_winds)){
      date = rep(current.date, each = length(current.idx))
      t.start = as.numeric(summary.current.RF$start[is.element(summary.current.RF$seg, current.idx)])
      t.stop = as.numeric(summary.current.RF$stop[is.element(summary.current.RF$seg, current.idx)])
      
      
      
      ###########   
      
      # WORKING WITH THE WRF-LES INPUTS

      t.start.object =  as.POSIXct(strptime(paste(current.date, floor(t.start), floor((t.start %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")
      t.start.idx = which.min(abs(Times - t.start.object)) - 15
      t.stop =  as.POSIXct(strptime(paste(current.date, floor(t.stop), floor((t.stop %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")
      t.stop.idx = which.min(abs(Times - t.stop)) + 15 
     # browser()
      if (t.stop.idx > length(Times)){
        t.stop.idx = length(Times)
      }

      rr.stack.crop <- rr.stack 
      rr.t <- subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5))
      
      rr.t.mean <- raster::calc( subset(rr.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = rms.cal)
      
      raster_center = c( (extent(rr.stack)[1] + extent(rr.stack)[2])/2, 
                         (extent(rr.stack)[3] + extent(rr.stack)[4])/2)
                      #c( (extent(e.inflow)[1] + extent(e.inflow)[2])/2, 
                      #   (extent(e.inflow)[3] + extent(e.inflow)[4])/2)
      cropped_center = c( (extent(source_raster)[1] + extent(source_raster)[2])/2, 
                          (extent(source_raster)[3] + extent(source_raster)[4])/2)
      inflow_correction = cropped_center - raster_center
      corrected_inflow = e.inflow
      corrected_inflow[1:2] = corrected_inflow[1:2] - inflow_correction[1]
      corrected_inflow[3:4] = corrected_inflow[3:4] - inflow_correction[2]
     # browser()
      e.inflow = corrected_inflow
      #browser()
      rr.inflow.stack.crop <- raster::crop(rr.stack, e.inflow)
      rr.inflow.t.mean <- raster::calc( subset(rr.inflow.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = rms.cal)
      
      psf.inflow.stack.crop <- raster::crop(psfc.stack, e.inflow)
      psd.inflow.t.mean <-  subset(psf.inflow.stack.crop,  seq(t.start.idx, t.stop.idx, 5))
      psd.inflow.t <- raster::calc( subset(psf.inflow.stack.crop,  seq(t.start.idx, t.stop.idx, 5)), fun = rms.cal)
      
      n_matched = sum(summary.current.source[, 'Tot.area'], na.rm = T)/median(raster::values(raster::area(rr.t.mean)) * 1e6)
      
      selected_clump= (raster::values( rr.t.mean ) > sort(raster::values(rr.t.mean) )[1])
      selected_clump[is.na(selected_clump)] = F
      #les.mass = sum( molef.to.kg ( (sort(raster::values(rr.t.mean) )[1:n_matched] ),  mean(raster::values(psd.inflow.t.mean), na.rm =  T)) ) * 0.111 * 0.111 # * raster::values(raster::area(rr.t.mean))[selected_clump] * 1e6, na.rm=T)
      les.mass.all = NULL
      
      for (idx_layer in seq(1, nlayers(rr.t))){
        x = subset(rr.t, names(rr.t)[idx_layer] )
        les.mass.all = c(les.mass.all, sum( molef.to.kg ( (sort(raster::values(x) , decreasing = T)),  mean(raster::values(psd.inflow.t), na.rm =  T))* 0.111 * 0.111, na.rm=T))
      }
      obs.mass =  sum(summary.current.source[, 'IME.kg'], na.rm = T)
      
      ratio.estimate.all = obs.mass/les.mass.all * 988
      
      ratio.sd =  sd(ratio.estimate.all)
      
      bt.mean = bootstrap::bootstrap(ratio.estimate.all, 1000, rms.cal)$thetastar
      ratio.estimate = quantile(bt.mean, 0.5) 
      ratio25 = quantile(bt.mean, 0.025) #qs.emisss[1]
      ratio975 = quantile(bt.mean, 0.975) #qs.emisss[2]
      
      bt.prc = bootstrap::bootstrap(sqrt(((ratio.estimate.all - ratio.estimate)/ratio.estimate)**2)*100, 1000, rms.cal)$thetastar
      #+++++++++++++++++++++++++++++++
    
      wind.stack.crop <- raster::crop(wind.stack, e.inflow)
      
      #wind.stack.crop <- raster::crop(wind.stack, e.current)
      #browser()
      wind.stack.crop.t <-  subset(wind.stack.crop,  seq(t.start.idx, t.stop.idx, 5))
      
      # IME 
      nsamples = ncell(subset(wind.stack.crop.t, names(wind.stack.crop.t)[sample(1:((t.stop.idx - t.start.idx)/5 + 1), 1, replace = T)] ))
      current.wind = rms.cal(array(values(wind.stack.crop.t))[sample(1:length(array(values(wind.stack.crop.t))), nsamples * ((t.stop.idx - t.start.idx)/5 + 1), replace = T)])
      current.winds = raster::values(wind.stack.crop.t) #cellStats(wind.stack.crop.t, "rms")
      current.ime =  sum(summary.current.source[, 5][summary.current.source[, 1] ==  current.segment.input])
      current.imes =  summary.current.source[, 5][summary.current.source[, 1] ==  current.segment.input]
      current.areas  = summary.current.source[, 7][summary.current.source[, 1] ==  current.segment.input & summary.current.source[, 5] > 0]
      current.area  = sum(summary.current.source[, 7][summary.current.source[, 1] ==  current.segment.input & summary.current.source[, 5] > 0])
     #browser()
      ime.emiss = emission.cal.IME(wind = current.wind, ime_input = current.ime, area_input = current.area)
      print("===================================")
      print("IME (kg/hr)")
      print(ime.emiss)
      
      # if(ime.emiss < 0){
      #   ime.emiss = NA #browser()
      # }
    } else{
      current.wind =  skip_winds
      
      current.ime =  sum(summary.current.source[, 5][summary.current.source[, 1] ==  current.segment.input])
      current.imes =  summary.current.source[, 5][summary.current.source[, 1] ==  current.segment.input]
      current.areas  = summary.current.source[, 7][summary.current.source[, 1] ==  current.segment.input & summary.current.source[, 5] > 0]
      current.area  = sum(summary.current.source[, 7][summary.current.source[, 1] ==  current.segment.input & summary.current.source[, 5] > 0])
      print("Calculating the emissions and quantiles")
      ime.emiss = emission.cal.IME(wind = current.wind, ime_input = current.ime, area_input = current.area)
      min25 = 0#quantile(bt.mean, 0.025) #qs.emisss[1]
      max975 = 0#quantile(bt.mean, 0.975) #qs.emisss[2]
      sd5 =  0 #sd(ime.emisss)
      print(ime.emiss)
      ratio.estimate = 0
      ratio.sd = 0
      ratio25 = 0 
      ratio975 = 0
      
    }
      # END SKIPPING WINDS
      
      print("===================================")
      mean.signal = mean(summary.current.source[, 'signal'])
      sd.noise = mean(summary.current.source[, 'noise'] )
      SNR = mean.signal/sd.noise
      
      min25 = 0
      max975 = 0
      sd5 = 0
      
      potential_output = list(current.segment.input, ime.emiss,  
                              current.area, current.ime, 
                              length(current.imes),  
                              ratio.estimate, 
                              ratio.sd, ratio25, 
                              ratio975, min25, 
                              max975, sd5, 0, 
                              0, current.wind, 
                              current_thres, SNR)
      
      if (sum(is.na(potential_output)) == 0){
        df_row = data.frame(potential_output)
        names(df_row)<- names(current_row_df_summary)
        current_row_df_summary <- rbind(current_row_df_summary, df_row)
        
      }else (print('skip this clump'))
      
     # browser()
      
      print(paste("Done with this execution: ", current.research.flight, source.name, as.character(current.segment.input)))
      print("==============================================================")
      current_row_df_summary %>% mutate_at(vars(IME.EMISS, IME.AREA, IME.IME, RATIO.EMISS, RATIO.SD, RATIO.25, RATIO.975, IME.25, IME.975, SD, X.WIND), function(x)(round(x, 2))) -> current_row_df_summary
      rownames(current_row_df_summary) <- seq(dim(current_row_df_summary)[1])
      
     # browser()
      if(exists("t.start")){
      current_row_df_summary$formatted.time  <- as.POSIXct(strptime(paste(current.date, floor(t.start), floor((t.start %% 1)  * 60)), "%Y-%m-%d %H %M"), tz = "UTC")
      } else{current_row_df_summary$formatted.time = NA}
      gc()
    }
    else{current_row_df_summary = NULL
    les.mass.all = NULL
    current.winds = NULL}
    
    df_summary_combined = rbind(df_summary_combined, current_row_df_summary)
  # browser()
    ######################## ABOVE == ORIGINAL mIME ########################
    ######################## BELOW == RESAMPLING mIME ########################
  # LOOPING OVER THE N_i SAMPLES   
  # START RESAMPLING   (MCMC)

    if(is.null(resample_estimates) && !skipping){
    print("START RESAMPLING")
    spatial_points = as(source_above_thres>0, "SpatialPoints")
    crs(spatial_points) <- crs(source_above_thres)
    buffer_object = gBuffer(spatial_points, width = 0.0001)
   # browser()
    temp_source_above_thres_radius = mask(source_above_thres, buffer_object)
    
    
    source_array = values(mask(source_raster, buffer_object))
    source_array = source_array[!is.na(source_array)] 
    inflow_array = values(source_inflow)
    resample.ime.emisss = c()
    resample_kgs = c()
    

  for (idx_resample_thre in seq(N_i)) { 
   # print(idx_resample_thre)
    selected_clumps = NULL
   
    resample_source = raster::calc(source_raster, fun = function(x) sample(source_array, 1, replace = T))
    resample_source = mask(resample_source, buffer_object)

    #plot(resample_source)
    # plot(resample_source, col = heat.275, main = "resample source")
    #browser()
    resample_inflow = sample(inflow_array, size = length(inflow_array), replace = T)
    
#    plot(resample_inflow, col = heat.275, main = "resample inflow")
    # test_raster = source_inflow
    # test_raster$xch4 = resample_inflow
    # plot(test_raster, col = heat.275)
    # browser()
    
    resample_inflow_median = median( resample_inflow,  na.rm=T ) 
    resample_inflow_sd =  sd(resample_inflow, na.rm = T)
    resample_threshold = resample_inflow_median + resample_inflow_sd * current_thres
    
    
    resample_source_above_thres = resample_source
    raster::values(resample_source_above_thres)[raster::values(resample_source_above_thres) < resample_threshold] = NA
    # print("LOOK AT THIS")
    # plot(resample_source_above_thres, col = heat.275, main = "resample above thre")
    # browser()
    
    ### save the median, sd, and threshold from each iteration
      my.median = c(my.median, inflow_median) 
      my.sd = c(my.sd, inflow_sd) 
      my.thr = c(my.thr,  resample_threshold)
    ### save the median, sd, and threshold from each iteration
    resample_kg = sum( molef.to.kg ( (raster::values(resample_source_above_thres) - resample_inflow_median ), 
                            entire_pressure_input) * raster::values(raster::area( resample_source))[values(resample_source) > 0] * 1e6, na.rm=T)
    if(is.null(skip_winds)){    
      current.wind = rms.cal(array(values(wind.stack.crop.t))[sample(1:length(array(values(wind.stack.crop.t))), nsamples * ((t.stop.idx - t.start.idx)/5 + 1), replace = T)])
      }else{current.wind = current.wind}
  # browser()
   
   new_area = sum(raster::values(raster::area(resample_source_above_thres))[values(resample_source_above_thres) > 0], na.rm=T) * 1e6
   resample.ime.emiss = emission.cal.IME(wind = current.wind, ime_input = resample_kg, area_input = new_area, coef = runif(1, 0.9, 1.1))
   
    #resample.ime.emiss = emission.cal.IME(wind = current.wind, ime_input = resample_kg, area_input = current.area, coef = runif(1, 0.9, 1.1))
    
    resample_kgs = c(resample_kgs, resample_kg)
    resample.ime.emisss = c(resample.ime.emisss, resample.ime.emiss)
   ########## NEW LOOP OF MCMC ################## 
    print(idx_resample_thre)
  } 
    #browser()
  #print(resample.ime.emisss)
    df_summary_output = current_row_df_summary
    lower_CI = function(x) quantile(x, 0.025, na.rm = T)
    upper_CI = function(x) quantile(x, 0.975, na.rm = T)
    median_CI = function(x) quantile(x, 0.5, na.rm = T)
    df_summary_output$IME.25  = lower_CI(resample.ime.emisss) #lower_CI(bootstrap::bootstrap(resample.ime.emisss, N_i, median)$thetastar)
    df_summary_output$IME.975 = upper_CI(resample.ime.emisss) #upper_CI(bootstrap::bootstrap(resample.ime.emisss, N_i, median)$thetastar)
    df_summary_output$IME.EMISS = median_CI(resample.ime.emisss) #median_CI(bootstrap::bootstrap(resample.ime.emisss, N_i, median)$thetastar)
    df_summary_output$IME.IME = median(resample_kgs, na.rm = T) 
#browser()
    
 if(save_raw){
   raw_row = data.frame(RF = current.research.flight, seg = current.segment.input, 
                        formatted.time = current_row_df_summary$formatted.time, rawbt = resample.ime.emisss)
 }else{
   raw_row = NULL
 }

    current.date = date(current_row_df_summary$formatted.time)
    YYYYMMDD = str_replace_all(current.date, "-", "")
    seconds_since_midnight =  hour(current_row_df_summary$formatted.time) * 60 * 60 + minute(current_row_df_summary$formatted.time) * 60
    export_name = paste(YYYYMMDD, seconds_since_midnight, format(raw_maryann_data$x_plume[current.segment.input], nsmall = 6), format(raw_maryann_data$y_plume[current.segment.input], nsmall = 6), "RF06_mIME.nc", sep = "_")
    print(export_name)
   # browser()
     if(export_data){
          
          export_nc = stack(source_raster, source_above_thres)
          names(export_nc)  = c("original", "cropped_bg_removed")
          current_output_nc = export_name
          current_output_nc_og =  paste("original_", export_name, sep = "")
          current_output_nc_cropped =  paste("cropped_", export_name, sep = "")
          current_output_mask =  paste("masked_", export_name, sep = "")
          #browser()
          writeRaster(export_nc, filename=file.path(out.netcdf, current_output_nc), format="CDF", overwrite=TRUE)   
          names(source_raster) = "xch4"
          writeRaster(source_raster, filename=file.path(out.netcdf, current_output_nc_og), format="CDF", overwrite=TRUE)
          
          names(source_above_thres) = "cropped_xch4"
          writeRaster(source_above_thres, filename=file.path(out.netcdf, current_output_nc_cropped), format="CDF", overwrite=TRUE)   
          
          masked_raster = source_above_thres > 0
          names(masked_raster) = "mask"
          writeRaster(masked_raster , filename=file.path(out.netcdf, current_output_mask), format="CDF", overwrite=TRUE)   
         # browser()
          KML(source_raster, filename=file.path(out.netcdf, current_output_nc_og), col=rev(terrain.colors(255)), 
              colNA=NA, maxpixels=100000, blur=1, zip='', overwrite=TRUE)
          KML(source_above_thres, filename=file.path(out.netcdf, current_output_nc_cropped), col=rev(terrain.colors(255)), 
              colNA=NA, maxpixels=100000, blur=1, zip='', overwrite=TRUE)
     } else{export_name = NULL}
    
    # !!!!!  MIGHT WANT TO EXPORT df_summary_combined  !!!!!
    # So we can refer to the original bootstrapped 1000 samples 
#browser()
  return(list(df_summary_output[1, ],  list(les.mass.all, as.vector(current.winds), export_name, raw_row) ))}
    else{
      return(list(df_summary_combined[1, ],  cbind(NULL, NULL, export_name, raw_row) ))
    }
}

