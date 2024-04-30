# This script takes in list of lat/lons of plumes, and finds the upwind end, plume length, 
# and plume-derived wind angle. Runs DI growing box on each plume.

setwd("~/Dropbox/MethaneAIR")

library(ncdf4)
library(fields)
library(raster)
library(class)
library(geosphere)
library(DescTools)
library(smatr)


# Inputs
RF="06"
date = "20210806"       # flight date format YYYYMMDD  
hour=16        
# lat/lons of plumes
#xy.plume1=rbind(c(-104.7616,39.851),c(-105.013, 40.0228),c(-104.832,40.5854),c(-104.5,40.367),c(-104.783,40.325)) #tower; front range; north weld; five rivers Kerner; five rivers Gilcrest
aa=read.csv("L4 plume flags/RF06_po444b_QAQC_plumes.csv",sep=",")
xy.plume1=as.matrix(aa[,c("Lon","Lat")])

file.typ="MRS" # DPP or MRS
l3.file="L3filtered/MSAT_RF06_L3_filt.85.444b.nc"

# Adjustable parameters
sig.sel=1.5
nmin=120 # typical ~40 for 30m; 100 for 10m filtered

########
# force script to use particular locations, plume lengths, or wind rotation for certain flights:
# note, if you want to change plume location by hand, change lat/lon in xy.plume1 by hand
plume.loc=rep(2,dim(xy.plume1)[1]) # 1 = use input plume location; 2 = use plume-derived plume location; 3=use whichever is upwind between input and plume-derived
plume.length=rep(NA,dim(xy.plume1)[1]) # NA = use plume-derived plume length; # = use # for plume length
rotate_force=rep(-2,dim(xy.plume1)[1]) # -2 = use script below to rotate wind or not; -1= never rotate wind (use HRRR); # (0-360): use # for wind angle

if(RF=="06"){
  plume.loc[c(8,9,10,11)]=1
}


#################################

source("emissions_DI_growing_box2.R")
load("heat.275.RData")

if(file.typ=="DPP"){
  xch4=raster( l3.file ,varname="xch4")
  num_samples=raster(my_mosaic_local_files[lf], varname = "num_samples")
  values(xch4)[values(xch4)>1e10 | values(num_samples) < 0.2  ]=NA
  values(xch4)=values(xch4)*1E-9
} else{
  xch4=raster(l3.file) # L3 file in raster format 
  values(xch4)=values(xch4)*1E-9
  
}

dat1a=as.data.frame(xch4, xy = TRUE)
met=load.met(date, hour, dat1a,RF) 

flux=NULL

for(i in 1:dim(xy.plume1)[1]){ 
  
  xy.plume2=xy.plume1[i,]
  xch4a=crop(xch4,c(xy.plume1[i,1]-.02,xy.plume1[i,1]+.02,xy.plume1[i,2]-.02,xy.plume1[i,2]+.02))
  dat1=as.data.frame(xch4a, xy = TRUE)
  

    
  # thresholds to find xch4 plumes
  my.trim.xch4 = median(values(xch4a),na.rm=T) + sig.sel * sd(values(xch4a),na.rm=T) # threshold 1 used for which plumes to keep
  my.trim.xch4.2 = median(values(xch4a),na.rm=T) + sig.sel*1.3 * sd(values(xch4a),na.rm=T)  # threshold 2 used for deciding whether to rotate wind to match plume direction
  
  # new xch4 rasters including only values above threshold
  xch4a.sel = xch4a
  values( xch4a.sel ) [ values( xch4a) < my.trim.xch4 ] = NA  
  my.clumps.xch4 = clump( xch4a.sel,directions=8)   # find clumps
  
  # repeat with higher threshold: new xch4 rasters including only values above higher threshold
  xch4.plume.sel.2 = xch4a
  values( xch4.plume.sel.2 ) [ values( xch4a) < my.trim.xch4.2 ] = NA  
  my.clumps.xch4.2 = clump( xch4.plume.sel.2,directions=8)    # find clumps
  
  
  # intialize raster to hold only selected xch4 clumps (above threshold number of pixels)
  xch4.final = xch4a.sel
  values(xch4.final) = NA    
  xch4.final.2 = xch4a.sel
  values(xch4.final.2) = NA       
  
  
  nclumps.xch4 = unique( values( my.clumps.xch4 ))
  nclumps.xch4 = nclumps.xch4[ is.finite( nclumps.xch4) ]
  vv.xch4 = values(my.clumps.xch4)
  pp.xch4 = values( xch4a.sel )
  xy.xch4 = coordinates ( xch4a )
    
  # Only keep large clumps (#pixels > nmin.xch4).  Create xch4.final holding all clumps above threshold and above min size (all else=NA)
  # to keep a clump, part of that clump must be within 1 km of the center of the DI plume (this is relevant for fitting a line to the plumes later)
  for ( m in 1:length(nclumps.xch4)) {
    l.selx = ( vv.xch4 == m )
    l.selx[ is.na(l.selx) ] = F
    if (sum ( l.selx ) < nmin ) next # only include clumps with more cells than nmin.xch4
    dist1=min(distm(xy.xch4[l.selx,],xy.plume2)) # distance from DI clump center to closest part of xch4 clump
    if(dist1<1000){values(xch4.final)[l.selx ] = pp.xch4[ l.selx ]} # keep xch4 values of selected clump
  }
  
  
  
  # Repeat clumping with higher threshold
  nclumps.xch4.2 = unique( values( my.clumps.xch4.2 ))
  nclumps.xch4.2 = nclumps.xch4.2[ is.finite( nclumps.xch4.2) ]
  vv.xch4.2 = values(my.clumps.xch4.2)
  pp.xch4.2 = values( xch4.plume.sel.2 )
  xy.xch4.2 = coordinates ( xch4.plume.sel.2 )
  
  # create xch4.final with xch4 in clumps above threshold number of pixels (nmin.xch4), all greater than my.trim.value.xch4, and other cells = NA
  for ( m in 1:length(nclumps.xch4.2)) {
    l.selx = ( vv.xch4.2 == m )
    l.selx[ is.na(l.selx) ] = F
    if (sum ( l.selx ) < nmin ) next # only include clumps with more cells than nmin
    values(xch4.final.2)[l.selx ] = pp.xch4.2[ l.selx ]
  }
  # to label a strong plume, there must be part of the higher threxhold xch4 plume within .01o (1km) of center of DI plume 
  xch4.final2.2=crop(xch4.final.2,c(xy.plume2[1]-.01,xy.plume2[1]+.01,xy.plume2[2]-.01,xy.plume2[2]+.01)) # crop in closer to original flux plume
  
    
if(sum(!is.na(values(xch4.final)))>0){
    ####################
    # Find center of mass/moment of inertia of XCH4 plume
    cosfac= abs( cos( mean(xy.xch4[,2]) * pi/180) )   # cosine of the mean latitude:  longitude contraction factor for equidistant
    
    ## xy coordinates of masses, in equal distance metric
    X = xy.xch4[, 1] * cosfac
    X0.xch4 = sum ( (values(xch4.final)* X )[ !is.na(values(xch4.final)) ] ) / sum( values(xch4.final)[!is.na(values(xch4.final))] ) / cosfac  ## loc of COM, referred back to lon/lat
    
    Y = xy.xch4[, 2] 
    Y0.xch4 = sum ( (values(xch4.final)* Y )[ !is.na(values(xch4.final)) ] ) / sum( values(xch4.final)[ !is.na(values(xch4.final)) ] ) 
    
    ## define the moments of inertia
    Iyy =   sum(  values( xch4.final)[ !is.na(values(xch4.final)) ] *  (  X[!is.na(values(xch4.final))]  - X0.xch4* cosfac  )  ^2 )
    Ixx =   sum(  values( xch4.final)[ !is.na(values(xch4.final)) ] *  (  Y[ !is.na(values(xch4.final))]  - Y0.xch4  )^2 )
    Ixy = - sum(values( xch4.final)[ !is.na(values(xch4.final)) ]  * (  X[!is.na(values(xch4.final))]  - X0.xch4* cosfac  )*(  Y[ !is.na(values(xch4.final))]  - Y0.xch4  ) )
    
    ## --------------  find the principal moments and locate the axes with "angle"
    Mii = eigen(matrix( c(Ixx, Ixy, Ixy, Iyy),ncol=2,nrow=2))
    ## The eigenvvectors give the directions of the principal axes, and the 
    #         eigen­ values give the moments of inertia with respect to each of these axis.
    ## this is a planar figure, so there are only two that matter
    i.p = 2 - as.numeric( Mii$values[2] > Mii$values[1])  # the long axis is the smaller moment
    if( i.p == 2 ) i.p.1 = 1  else i.p.1 = 2
    angle= 180* atan2( Mii$vectors[ 2, i.p], Mii$vectors[1, i.p])/ pi
    
    Delta.x =  diff( range(  X [!is.na(values(xch4.final))] )/ cosfac )
    Delta.y =  Delta.x *tan(angle* pi/180)
    
    my.slope = Delta.y/Delta.x
    my.angle = atan2(Delta.x,Delta.y)*180/pi # wind angle based on plume MOI
    #  equivalent ellipse
    Ix1 = Mii$values[ i.p ]  
    Iy1 = Mii$values[ i.p.1 ]
    eccentricity = sqrt( 1 - Ix1/Iy1 )
    
    ########################
    
    
    # HRRR wind data at location of plume
    met.nn1=knn(train=met$met.latlon,test=c(X0.xch4,Y0.xch4),cl=factor(1:nrow(met$met.latlon))) # location of plume (X0,Y0) mapped to HRRR grid
    met.pt=met.square(met,met.nn1,c(X0.xch4,Y0.xch4))   # HRRR data at X0,Y0
    uwind=met.pt$uwind2 # U wind data at X0,Y0
    vwind=met.pt$vwind2 # V wind data at X0,Y0
    surf_ht=met.pt$surf_ht # surface height at X0,Y0
    
    # Use HRRR wind direction to find upwind/downwind ends of plume axis
    # if HRRR and plume-derived wind directions differ by >90o, then adjust plume direction by 180o
    my.angle.hrrr = atan2(uwind,vwind)*180/pi
    angl.diff=acos( sum(c(Delta.y,Delta.x)*c(vwind,uwind)) / ( sqrt(sum(c(Delta.y,Delta.x) * c(Delta.y,Delta.x))) * sqrt(sum(c(vwind,uwind) * c(vwind,uwind))) ) )*180/pi  #anglular distance between MOI angle and HRRR angle
    my.angle1=my.angle  # wind angle based on plume MOI
    if(angl.diff>90){ # if HRRR and MOI wind are in different directions, then MOI is in the wrong direction - swap by 180
      if(my.angle>0){my.angle1=my.angle-180}
      if(my.angle<=0){my.angle1=my.angle+180}
    }
    
    # find length of plume
    dat.1=as.data.frame(xch4.final, xy = TRUE) # convert XCH4 plume raster to data frame (values for plume, NA elsewhere)
    dat.2=data.frame(x=dat.1[!is.na(dat.1[,3]),1],y=dat.1[!is.na(dat.1[,3]),2]) # keep lat/lon for only non-NA XCH4 values (within plume) 
    
    Xl = seq(min(dat.2$x),max(dat.2$x),by=1E-6)
    Yl = my.slope*(Xl-X0.xch4)+Y0.xch4
    Xl=Xl[Yl<=max(dat.2$y)&Yl>=min(dat.2$y)]
    Yl=Yl[Yl<=max(dat.2$y)&Yl>=min(dat.2$y)]
    plume_l=distm(c(Xl[1],Yl[1]),c(Xl[length(Xl)],Yl[length(Yl)])) # plume length: distance from upwind to downwind edge of plume
    
       
  
  wspd=sqrt(uwind^2+vwind^2)
  # Rotate wind to match long axis of plume IF eccentricity>0.87 AND XCH4 is above higher threshold
  if(eccentricity>=.87 & sum(values(xch4.final2.2),na.rm=T)>0){
    uwind=wspd*sin(my.angle1*pi/180)
    vwind=wspd*cos(my.angle1*pi/180)
  }
  
  dat.1=as.data.frame(xch4.final, xy = TRUE) # convert XCH4 plume raster to data frame (values for plume, NA elsewhere)
  dat.2=data.frame(x=dat.1[!is.na(dat.1[,3]),1],y=dat.1[!is.na(dat.1[,3]),2]) # keep lat/lon for only non-NA XCH4 values (within plume) 
  
  # find best estimate of plume origin
  # plume lat/lon = whichever is farther upwind: either either input lat/lon or downwind end of XCH4 clump
  if(abs(uwind)>abs(vwind)){
    if(uwind>0){
      source3=c(min(dat.2$x),my.slope*(min(dat.2$x)-X0.xch4)+Y0.xch4) # source3: XCH4 plume upwind end
      if(source3[2]>max(dat.2$y)){source3[2]=max(dat.2$y)}
      if(source3[2]<min(dat.2$y)){source3[2]=min(dat.2$y)}
      source2=source3 # source2: plume origin lat/lon
      if(source3[1]>xy.plume2[1]){source2=c(xy.plume2[1],xy.plume2[2])}
      downwind3=c(max(dat.2$x),source2[2])  # downwind3: XCH4 downwind end (for calculating plume length) (as of now, calculate length along X or Y depending on wind direction, not diagonal (if we switch to rotated boxes, would want to use diagonal))
    } else if(uwind<=0){
      source3=c(max(dat.2$x),my.slope*(max(dat.2$x)-X0.xch4)+Y0.xch4) # source3: XCH4 plume upwind end
      if(source3[2]>max(dat.2$y)){source3[2]=max(dat.2$y)}
      if(source3[2]<min(dat.2$y)){source3[2]=min(dat.2$y)}
      source2=source3 # source2: plume origin lat/lon
      if(source3[1]<xy.plume2[1]){source2=c(xy.plume2[1],xy.plume2[2])}
      downwind3=c(min(dat.2$x),source2[2]) # downwind3: XCH4 downwind end (for calculating plume length) (as of now, calculate length along X or Y depending on wind direction, not diagonal (if we switch to rotated boxes, would want to use diagonal))
    }
  } else if(abs(uwind)<=abs(vwind)){
    if(vwind>0){
      source3=c((min(dat.2$y)-Y0.xch4)/my.slope+X0.xch4,min(dat.2$y))  # source3: XCH4 plume upwind end
      if(source3[1]>max(dat.2$x)){source3[1]=max(dat.2$x)}
      if(source3[1]<min(dat.2$x)){source3[1]=min(dat.2$x)}
      source2=source3 # source2: plume origin lat/lon
      if(source3[2]>xy.plume2[2]){source2=c(xy.plume2[1],xy.plume2[2])}
      downwind3=c(source2[1],max(dat.2$y)) # downwind3: XCH4 downwind end (for calculating plume length) (as of now, calculate length along X or Y depending on wind direction, not diagonal (if we switch to rotated boxes, would want to use diagonal))
    } else if(vwind<=0){
      source3=c((max(dat.2$y)-Y0.xch4)/my.slope+X0.xch4,max(dat.2$y))  # source3: XCH4 plume upwind end
      if(source3[1]>max(dat.2$x)){source3[1]=max(dat.2$x)}
      if(source3[1]<min(dat.2$x)){source3[1]=min(dat.2$x)}
      source2=source3 # source2: plume origin lat/lon
      if(source3[2]<xy.plume2[2]){source2=c(xy.plume2[1],xy.plume2[2])}
      downwind3=c(source2[1],min(dat.2$y)) # downwind3: XCH4 downwind end (for calculating plume length) (as of now, calculate length along X or Y depending on wind direction, not diagonal (if we switch to rotated boxes, would want to use diagonal))
    }
  }
  
  # keep output data for plume
  centXch4=c(source3[1],source3[2])
  centAll=c(source2[1],source2[2])
  rotate_wind=eccentricity>=.87 & sum(values(xch4.final2.2),na.rm=T)>0 # if T, rotate wind; if F, use HRRR wind angle
  ########

  par(mfrow=c(1,2),mar=c(5,4.5,2,1))
  plot(xch4a*1E9,main=i,col=heat.275)
  plot(xch4.final*1E9,main=paste(i,round(my.angle,2)),col=heat.275)
  if(sum(!is.na(values(xch4.final)))>0){
    lines( Xl, Yl,col=1)}
  points(xy.plume2[1],xy.plume2[2],col=2)
  points(centXch4[1],centXch4[2],col=4)
  points(centAll[1],centAll[2],pch=3)
  legend('topright',c("input","plume","upwind"),col=c(2,4,1),pch=c(1,1,3))
  
} else{
  my.angle1=NA
  my.angle="HRRR"
  plume_l=800
  centXch4=xy.plume2
  centAll=xy.plume2
  rotate_wind=0
}
  
  if(plume.loc[i]==1){xy.final=xy.plume2}
  if(plume.loc[i]==2){xy.final=centXch4}
  if(plume.loc[i]==3){xy.final=centAll}
  
  if(is.na(plume.length[i])){plume_len=as.numeric(plume_l)
  }else{
    plume_len=plume.length[i]
  }
  if(rotate_force[i]==(-1)){rotate_wind=0}
  if(rotate_force[i]>=0){
    rotate_wind=1
    my.angle1=rotate_force[i]}
  
  wdir.obs="HRRR"
  if(rotate_wind==1){wdir.obs=my.angle1}
  if(wdir.obs<0){wdir.obs=wdir.obs+360}
  

  
    f1=emissions_DI_growing_box2(
      RF=RF,                                # Flight number
      file.typ="SCW",                          # Options: "CCM" for files by Chris Miller of form: regridded_diag_0p85pixwidth_vida-pass1.nc; "SCW" for files by Steve Wofsy of form: "RF06_mosaic.nc"; "MRS" for files by Maryann Sargent of form: "xch4_gaussfilt.RF07.segment.9.test.RData" 
      wdir.obs= wdir.obs,                        # if ="MOI", loads file with wind directions in scene calculated using Moment of Inertia method (moments.fcn3.R); if="HRRR", wind direction taken from HRRR; if=number, single wind direction for whole scene (wind from south to north=0o, wind from west to east=90o, etc.)
      my.raster=xch4,                # paste("~/Dropbox/MethaneAIR/RF06 SEGMENT_RASTERS/xch4.RF06.segment.",seg,".nc",sep="") # "~/Dropbox/MethaneAIR/rf08_composite_files/xch4.scr.agg10_rf08_uinta_mosaic_trimmed.nc"  ##"RF06_mosaic.nc"          # gridded "xch4" scene file     ex:"xch4_gaussfilt_RF07_5x1_mosaic.nc"  #  #"~/Dropbox/MethaneAIR/rf08_composite_files/xch4.scr.agg10_rf08_uinta_mosaic_trimmed.nc" #   #"xch4_bias_corr_v2_rf08_uinta_mosaic_trimmed.nc"  #"xch4_gaussfilt_RF07_5x1_mosaic.nc",  #"RF06_mosaic.nc"
      dat1=as.data.frame(xch4, xy = TRUE),
      date = date,                      # flight date format YYYYMMDD  #RF06 806; RF07 809
      xy.plume= xy.final,   # c(-104.756,39.8505), #c(-111.78591,32.824),  #c(-103.47244,   31.93576)  #c(-109.46,39.93)  #c(-109.37,40.15)  #c(-103.6, 32.1)  #centAll[25,]           #c(-103.6,32.1)                # source location    c(-109.78,39.95)  #c(-109.37,40.15)   # #c(-109.78,40)                  # lon/lat corners of box for flux calculations
      plume_l=plume_len,   # length of plume in m
      hour=hour,                                # hour (UTC) closest to measurement time.  RF06: seg1-3 16h; seg4-10 17h; 11-18 18h.  RF07: seg0-1 16h; seg2-4 17h; seg5 18h
      plot.on=T,                             # if TRUE, plot xch4 scene and tiles
      name.out="old",
      plume_num=i,
      my.box=c((xy.final[1]-(plume_l*9E-6+.035)),(xy.final[1]+(plume_l*9E-6+.035)),(xy.final[2]-(plume_l*9E-6+.035)),(xy.final[2]+(plume_l*9E-6+.035)))   #      # lon/lat corners of box for to bound xch4 map
    )
    
    flux=rbind(flux,c(i,xy.plume2,centXch4, centAll,plume.loc[i],f1$DI.mn, rotate_wind,my.angle1,plume_l))
  
}

flux2=flux[,-(12:16)]
colnames(flux2)=c("plume #","lat input","lon input","lat plume","lon plume","lat upwind","lon upwind","xy used","DI flux","flux SD","flux LCI","flux HCI",colnames(flux)[17:21],"rotate plume?","plume angle","plume length")


