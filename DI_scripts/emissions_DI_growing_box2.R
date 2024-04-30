# Function to calculate flux divergence for MethaneAir point source or regional sources (~1 km - 30 km scale) using Divergence Integral (Gauss theorum) method


emissions_DI_growing_box2= function(
    
  ################
  # Input parameters
  RF="08",                                # Flight number
  file.typ="SCW",                          # Options: "CCM" for files by Chris Miller of form: regridded_diag_0p85pixwidth_vida-pass1.nc; "SCW" for files by Steve Wofsy of form: "RF06_mosaic.nc"; "MRS" for files by Maryann Sargent of form: "xch4_gaussfilt.RF07.segment.9.test.RData" 
  wdir.obs= "HRRR",                        # if ="MOI", loads file with wind directions in scene calculated using Moment of Inertia method (moments.fcn3.R); if="HRRR", wind direction taken from HRRR; if=number, single wind direction for whole scene (wind from south to north=0o, wind from west to east=90o, etc.)
  my.raster=dum$my.raster,                # paste("~/Dropbox/MethaneAIR/RF06 SEGMENT_RASTERS/xch4.RF06.segment.",seg,".nc",sep="") # "~/Dropbox/MethaneAIR/rf08_composite_files/xch4.scr.agg10_rf08_uinta_mosaic_trimmed.nc"  ##"RF06_mosaic.nc"          # gridded "xch4" scene file     ex:"xch4_gaussfilt_RF07_5x1_mosaic.nc"  #  #"~/Dropbox/MethaneAIR/rf08_composite_files/xch4.scr.agg10_rf08_uinta_mosaic_trimmed.nc" #   #"xch4_bias_corr_v2_rf08_uinta_mosaic_trimmed.nc"  #"xch4_gaussfilt_RF07_5x1_mosaic.nc",  #"RF06_mosaic.nc"
  dat1=dum$dat1,
  date = "20210811",                      # flight date format YYYYMMDD  #RF06 806; RF07 809
  xy.plume= plume.info1[i,c("x_plume","y_plume")],  #c(-103.47244,   31.93576)  #c(-109.46,39.93)  #c(-109.37,40.15)  #c(-103.6, 32.1)  #centAll[25,]           #c(-103.6,32.1)                # source location    c(-109.78,39.95)  #c(-109.37,40.15)   # #c(-109.78,40)                  # lon/lat corners of box for flux calculations
  plume_l=plume.info1[i,"plume_length"],   # length of plume in m
  hour=18,                                # hour (UTC) closest to measurement time.  RF06: seg1-3 16h; seg4-10 17h; 11-18 18h.  RF07: seg0-1 16h; seg2-4 17h; seg5 18h
  plot.on=T,                             # if TRUE, plot xch4 scene and tiles
  name.out="1.35",
  plume_num=1,
  my.box=c((xy.plume[1]-.05),(xy.plume[1]+.05),(xy.plume[2]-.05),(xy.plume[2]+.05))   #      # lon/lat corners of box for to bound xch4 map
  ################
){

setwd("~/Dropbox/MethaneAIR")
  

  library(ncdf4)
  library(raster)
  library(geosphere)
  library(class)
  library(fields)
  library(DescTools)

  
  ################
  # calculate flux for growing box method
  
  load("heat.275.RData")
  
  dat1=dat1[dat1$x>=my.box[1]&dat1$x<=my.box[2]&dat1$y>=my.box[3]&dat1$y<=my.box[4],]
  
  if(wdir.obs<0){wdir.obs=wdir.obs+360}
  
  if(plot.on==T){
    my.raster2=my.raster
    values(my.raster2)=values(my.raster)*1E6
    #my.raster3=crop(my.raster2,extent(c(my.box[1]+.025,my.box[2]-.025,my.box[3]+.025,my.box[4]-.025)))
    my.raster3=crop(my.raster2,extent(c(my.box[1],my.box[2],my.box[3],my.box[4])))
    pdf(paste("RF",RF,"_DI_plume",plume_num,"_",round(xy.plume[2],2),"_",-1*round(xy.plume[1],2),".pdf",sep=""), width = 4.8, height = 2.6, pointsize=7.5,paper = "special", compress=F,colormodel = "rgb",font="Helvetica",onefile=F)
    par(mfrow=c(1,2),mar=c(5,4.5,2,3))
    plot(my.raster3,col=heat.275,xlim=c(my.box[1]+.025,my.box[2]-.025),ylim=c(my.box[3]+.025,my.box[4]-.025))
    #dev.off()
    
    
  }
  # load met file
  met=load.met(date, hour, dat1, RF)
  
  # find met data at lat/lon of plume location
  met.nn1=knn(train=met$met.latlon,test=xy.plume,cl=factor(1:nrow(met$met.latlon)))
  met.pt=met.square(met,met.nn1,cbind(xy.plume,xy.plume))
  wspd.pt=sqrt(met.pt$uwind2^2+met.pt$vwind2^2)
  
  wdir1=wdir.obs
  # calculate wind direction from uwind, vwind
  if(wdir.obs=="HRRR"){
    wdir1=180*atan(met.pt$uwind2/met.pt$vwind2)/pi
    if(wdir1<0 & met.pt$uwind2<0){wdir1=wdir1+360}
    if(wdir1<0 & met.pt$vwind2<0){wdir1=wdir1+180}
    if(wdir1>0 & met.pt$uwind2<0 & met.pt$vwind2<0){wdir1=wdir1+180}
  }
  print(wdir1)
  
  # convert XCH4 data from raster to matrix
  xx=xFromCol(my.raster)
  yy=yFromRow(my.raster)
  dat2 = raster::as.matrix(my.raster)
  
  # find xy.plume location in xch4 matrix
  xcent=which(abs(xx-xy.plume[1])==min(abs(xx-xy.plume[1])))
  ycent=which(abs(yy-xy.plume[2])==min(abs(yy-xy.plume[2])))
  # find distance in m of x and y xch4 grid
  ddx=as.numeric(distm(c(xx[xcent]-mean(diff(xx)),yy[ycent]),c(xx[xcent],yy[ycent])))
  ddy=as.numeric(distm(c(xx[xcent],yy[ycent]-mean(diff(yy))),c(xx[xcent],yy[ycent])))
  # calculate number of grid cells between xy.plume and edge of xch4 matrix in downwind direction
  if(abs(met.pt$uwind2)>=abs(met.pt$vwind2)){
    #max.i=min((dim(dat2)[2]-xcent),xcent)
    max.i=round(plume_l/ddx+200/ddx)
    dist=as.numeric(ddx)*seq(2,max.i,by=1) # create array of distances - at each distance will calculate DI for a box at that distance from the source 
  }
  if(abs(met.pt$uwind2)<abs(met.pt$vwind2)){
    #max.i=min((dim(dat2)[1]-ycent),ycent)
    max.i=round(plume_l/ddy+200/ddy)
    dist=as.numeric(ddy)*seq(2,max.i,by=1)
  }
  flux1=NA*dist
  
  ind=0
  inds=seq(2,max.i,by=1)
  
  # calculate DI for a box of length i grid cells around plume origin
  for(i in seq(2,max.i,by=1)){
    
    ind=ind+1
      
    sq= make.square2(wdir1, xx, yy, dat2, xcent, ycent,i)
      
    # Plot some boxes
    if(length(dist)>=29){
    if(is.element(ind,c(seq(25,(length(dist)-4),by=20)))){         # big box: c(10,25,50,80,110,140,170,200,230,260,290)   #seq(2,length(dist),by=5))){
      lines(sq$square[,"x"],sq$square[,"y"],col=2,lwd=.5)
    }}
    
    if(sum(is.na(sq$square[,3]))<(.2*dim(sq$square)[1])){ # only calculate flux if NAs are less than 20 of points on square
      
      #length between successive points along box in m
      d1=rep(NA,length(sq$square[,3]))
      for(i in 2:length(sq$square[,3])){
        d1[i-1]=distm(sq$square[(i-1),c(1,2)],sq$square[(i),c(1,2)])
      }
      d1[length(d1)]=distm(sq$square[1,c(1,2)],sq$square[dim(sq$square)[1],c(1,2)])
      
      # find met data along box
      met.nn=knn(train=met$met.latlon,test=sq$square[,c(1,2)],cl=factor(1:nrow(met$met.latlon)))
      if(wdir.obs=="MOI"){
        moi.nn=knn(train=moi.latlon,test=l1[,c(1,2)],cl=factor(1:nrow(moi.latlon))) # test=square[]
        moi.wdir=mean(moi.angl[moi.nn],na.rm=T)
      }
      
      # find met data along box
      met.sq=met.square(met,met.nn,sq$square)
      
      pblh=mean(as.vector(met$pblh0)[met.nn]) # platetary boundary layer height
      
      # combine relavent info about square to input to flux calculation
      square.info = list( 
        met.nn=met.nn,
        square=sq$square,
        d1=d1,
        u_fac=sq$u_fac,
        v_fac=sq$v_fac
      )
      
      # calculate DI flux for square
      flux1[ind] = flux.calc(square.info, met.sq) 
    } #if(sum(is.na(square[,spp]))<(.1*dim(square)[1])){
  }   
  
  # Decide over which distance values to average the flux
  # Criteria to use flux value:
    # 1) only start averaging after flux passes 50 kg/hr for the first time
    # 2) dist>100m
    # 3) dist<plume_l+d_step*2
    # 4) if the flux becomes negative and stays negative for the rest of the distances, exclude those negative values (keep negative values if the flux later becomes positive)
  flux.thresh=which(flux1>=50) #start averaging fluxes after the flux first passes a minimum threshold of 50 kg/hr
  flux2=flux1[dist>=100&dist>dist[flux.thresh[1]]&dist<=(plume_l+40)]
  dist2=dist[dist>=100&dist>dist[flux.thresh[1]]&dist<=(plume_l+40)]
  inds2=inds[dist>=100&dist>dist[flux.thresh[1]]&dist<=(plume_l+40)]
    
  negs=which(flux2<(median(flux2,na.rm=T)*.3))
  negsx=NULL
  for(n in negs){
    if(max(flux2[n:length(flux2)],na.rm=T)<(median(flux2,na.rm=T)*.3)){negsx=cbind(negsx,n)}
  }
  if(is.null(negsx)){negsx1=length(dist2)
  }else{negsx1=negsx[1]-1}
  
  if(is.na(median(flux2,na.rm=T))){
    flux3=NA
    dist3=NA
    inds3=inds
  } else{
    flux3=flux2[which(flux2>=(median(flux2,na.rm=T)*.3))[1]:negsx1]
    dist3=dist2[which(flux2>=(median(flux2,na.rm=T)*.3))[1]:negsx1]
    inds3=inds2[which(flux2>=(median(flux2,na.rm=T)*.3))[1]:negsx1]
  }
  
  ###########
  # 2nd threshold - keep more points on each side of plume:
  flux.thresha=which(flux1>=30) #start averaging fluxes after the flux first passes a minimum threshold of 50 kg/hr
  flux2a=flux1[dist>=70&dist>dist[flux.thresha[1]]&dist<=(plume_l+40)]
  dist2a=dist[dist>=70&dist>dist[flux.thresha[1]]&dist<=(plume_l+40)]
  inds2a=inds[dist>=70&dist>dist[flux.thresha[1]]&dist<=(plume_l+40)]
  
  negs=which(flux2a<0)
  negsx=NULL
  for(n in negs){
    if(max(flux2a[n:length(flux2a)],na.rm=T)<0){negsx=cbind(negsx,n)}
  }
  if(is.null(negsx)){negsx1=length(dist2a)
  }else{negsx1=negsx[1]-1}
  
  if(is.na(median(flux2a,na.rm=T))){
    flux3a=NA
    dist3a=NA
    inds3a=inds
  } else{
    flux3a=flux2a[which(flux2a>=(median(flux2a,na.rm=T)*.05))[1]:negsx1]
    dist3a=dist2a[which(flux2a>=(median(flux2a,na.rm=T)*.05))[1]:negsx1]
    inds3a=inds2a[which(flux2a>=(median(flux2a,na.rm=T)*.05))[1]:negsx1]
  }
  #######
  
  
  sq3a= make.square2(wdir1, xx, yy, dat2, xcent, ycent,inds3[1])
  sq3= make.square2(wdir1, xx, yy, dat2, xcent, ycent,inds3[length(inds3)])
  lines(sq3a$square[,c(1,2)],col=2,lwd=.5)
  lines(sq3$square[,c(1,2)],col=2,lwd=.5)
  sq4=c(min(sq3$square[,"x"]),max(sq3$square[,"x"]),min(sq3$square[,"y"]),max(sq3$square[,"y"]))
  
  mn.all=mean(flux1[dist>=100&dist<=(plume_l+40)],na.rm=T)
  mn=mean(flux3)
  sd1=sd(flux3)
  if(is.na(mn)){mn=mn.all
  sd1=sd(flux1[dist>=100&dist<=(plume_l+40)],na.rm=T)
  }
  mn.2=mean(flux3a)
  sd2=sd(flux3a)
  if(is.na(mn.2)){mn.2=mn.all
  sd2=sd(flux1[dist>=50&dist<=(plume_l+40)],na.rm=T)
  }
  DI.mn=c(mn,sd1,mn-sd1,mn+sd1,mn.2,sd2,mn.2-sd2,mn.2+sd2,mn.all,sq4,wspd.pt)
  names(DI.mn)=c("mn","sd","lci","hci","mn2","sd2","lci2","hci2","mn.all","box.xmin","box.xmax","box.ymin","box.ymax","wspd")
  mn3=mn
    
  if(plot.on==T){
    plot(dist,flux1,col=1,typ="l",lwd=1.5,xlab="",ylab="",main=paste(round(mn), " [",round(mn-sd1),", ",round(mn+sd1),"]",sep=""))
    #points(dist2,flux2,pch=16,cex=.5)
    #points(dist3a,flux3a,pch=16,cex=1,col=1)
    points(dist3,flux3,pch=16,cex=.8,col=1)
    title(xlab="Distance (m)", line=2.5, cex.lab=1.2)
    title(ylab="Flux (kg/hr)", line=2.2, cex.lab=1.2)
    #plot(dist/1000,flux1,col=1,typ="l",lwd=2,xlab="Dist (km)",ylab="Flux (kg/hr)")
    dev.off()
  }
  
  # key variables to return: "dist", "flux1" 
  #save(list=c("dist","flux1","DI.mn"),file=paste("RF",RF,"_",name.out,"_",xy.plume[1],"_",xy.plume[2],".RData",sep=""))    
  return( list(
    dist=dist3,
    flux1=flux3,
    DI.mn=DI.mn
  ))
  
}  
############################

                                      
  

## --------- load.source --------
# load gridded level 3 xch4 data; output as data frame
load.source = function(source_file, file.typ) {
  if(file.typ=="CCM"){
    source.dat=nc_open(source_file)
    lon0=ncvar_get(source.dat,"lon")
    lat0=ncvar_get(source.dat,"lat")
    xch4=ncvar_get(source.dat,"xch4_bias_corr")
    my.raster=raster(xch4[dim(xch4)[1]:1,])
    extent(my.raster) <- c(min(lon0),max(lon0),min(lat0),max(lat0))
    dat1=as.data.frame(my.raster, xy = TRUE)
  } else if(file.typ=="SCW"){
    my.raster=raster(source_file)
    dat1=as.data.frame(my.raster, xy = TRUE)
  } else if(file.typ=="MRS"){
    load(source_file)
    #my.raster=seg
    my.raster=NA
    dat1=dat3
  } else if(file.typ=="new"){
    load(source_file)
    #my.raster=seg
    my.raster=NA
  }
  
  
  colnames(dat1)=c(colnames(dat1)[1:2],"xch4")
  # crop scene
  
    
  return( list(
    dat1=dat1,
    my.raster=my.raster))
}
## ------ END function ---------


## ------------ load.met --------------
# load meteorological data; save relevant parameters in area of xch4 data
load.met = function(date, hour, dat1, RF) {

  date1=paste(substring(date,1,4),substring(date,5,6),substring(date,7,8),sep="")
  
  if(hour>=12&hour<=17){hh="_12-17"}
  if(hour>=18&hour<=23){hh="_18-23"}
  met0 = nc_open(paste(date1,hh,"_hrrr_",hour,".nc",sep=""))
  
  met.lon0=ncvar_get(met0,"LON")
  met.lat0=ncvar_get(met0,"LAT")
  xx=round(dim(met.lon0)[1]/2)
  yy=round(dim(met.lon0)[2]/2)
  
  met.lon=met.lon0[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"])]
  met.lat=met.lat0[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"])]
  pres=ncvar_get(met0,"PRES")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"]),]
  temp=ncvar_get(met0,"TEMP")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"]),]
  pres_surf0=ncvar_get(met0,"PRSS")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"])]
  pblh0=ncvar_get(met0,"PBLH")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"])]
  uwind=ncvar_get(met0,"UWND")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"]),]
  vwind=ncvar_get(met0,"VWND")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"]),]
  level_h0=ncvar_get(met0,"HGHT")[met.lon0[,yy]>=(min(dat1[,"x"]))&met.lon0[,yy]<=(max(dat1[,"x"])),met.lat0[xx,]>=min(dat1[,"y"])&met.lat0[xx,]<=max(dat1[,"y"]),]
  wspd=sqrt(uwind^2+vwind^2)
  
  pres_surf=as.vector(pres_surf0)
  
  met.latlon=cbind(as.vector(met.lon),as.vector(met.lat))
  
  return (list (
    met.latlon=met.latlon,
    level_h0=level_h0,
    pres=pres,
    uwind=uwind,
    vwind=vwind,
    wspd=wspd,
    pres_surf=pres_surf,
    temp=temp,
    pblh0=pblh0
  ))
} 
## ----------- END function -----------

##-------------make.square2------------

make.square2 = function(wdir1, xx, yy, dat2, xcent, ycent,i) {
  # box sides of length i in 3 directions;  side of length i/4 in upwind direction
  # wind direction - angle from N
  if(wdir1>=225 & wdir1<=315){
    # create box of side length i
    l1x=(xcent-i):ceiling(xcent+i/4)
    l1y=rep((ycent-i),length(l1x))
    l2y=(ycent-i):(ycent+i)
    l2x=rep(ceiling(xcent+i/4),length(l2y))
    l3x=ceiling(xcent+i/4):(xcent-i)
    l3y=rep((ycent+i),length(l3x))
    l4y=(ycent+i):(ycent-i)
    l4x=rep((xcent-i),length(l2y))
    # lat/lon coordinages of box
    box=cbind(xx[c(l1x[-1],l2x[-1],l3x[-1],l4x[-1])],yy[c(l1y[-1],l2y[-1],l3y[-1],l4y[-1])])
    # xch4 around box
    sq1=c(dat2[(ycent-i),l1x[-1]],dat2[l2y[-1],ceiling(xcent+i/4)],dat2[(ycent+i),l3x[-1]],dat2[l4y[-1],(xcent-i)])
    square=cbind(box,sq1)
    colnames(square)=c("x","y","xch4")
  }
  if(wdir1>=45 & wdir1<=135){ 
    # create box of side length i
    l1x=ceiling(xcent-i/4):(xcent+i)
    l1y=rep((ycent-i),length(l1x))
    l2y=(ycent-i):(ycent+i)
    l2x=rep((xcent+i),length(l2y))
    l3x=(xcent+i):ceiling(xcent-i/4)
    l3y=rep((ycent+i),length(l3x))
    l4y=(ycent+i):(ycent-i)
    l4x=rep(ceiling(xcent-i/4),length(l2y))
    # lat/lon coordinages of box
    box=cbind(xx[c(l1x[-1],l2x[-1],l3x[-1],l4x[-1])],yy[c(l1y[-1],l2y[-1],l3y[-1],l4y[-1])])
    # xch4 around box
    sq1=c(dat2[(ycent-i),l1x[-1]],dat2[l2y[-1],(xcent+i)],dat2[(ycent+i),l3x[-1]],dat2[l4y[-1],ceiling(xcent-i/4)])
    square=cbind(box,sq1)
    colnames(square)=c("x","y","xch4")
  }
  if(wdir1>315 | wdir1<45){
    # create box of side length i
    l1x=(xcent-i):(xcent+i)
    l1y=rep((ycent-i),length(l1x))
    l2y=(ycent-i):ceiling(ycent+i/4)
    l2x=rep((xcent+i),length(l2y))
    l3x=(xcent+i):(xcent-i)
    l3y=rep(ceiling(ycent+i/4),length(l3x))
    l4y=ceiling(ycent+i/4):(ycent-i)
    l4x=rep((xcent-i),length(l2y))
    # lat/lon coordinages of box
    box=cbind(xx[c(l1x[-1],l2x[-1],l3x[-1],l4x[-1])],yy[c(l1y[-1],l2y[-1],l3y[-1],l4y[-1])])
    # xch4 around box
    sq1=c(dat2[(ycent-i),l1x[-1]],dat2[l2y[-1],(xcent+i)],dat2[ceiling(ycent+i/4),l3x[-1]],dat2[l4y[-1],(xcent-i)])
    #if(seg==9&ind>=50){
    #  l1x=(xcent-50):(xcent+i)
    #  l1y=rep((ycent-i),length(l1x))
    #  l2y=(ycent-i):ceiling(ycent+i/4)
    #  l2x=rep((xcent+i),length(l2y))
    # l3x=(xcent+i):(xcent-50)
    #  l3y=rep(ceiling(ycent+i/4),length(l3x))
    #  l4y=ceiling(ycent+i/4):(ycent-i)
    #  l4x=rep((xcent-50),length(l2y))
    #  # lat/lon coordinages of box
    #  box=cbind(xx[c(l1x[-1],l2x[-1],l3x[-1],l4x[-1])],yy[c(l1y[-1],l2y[-1],l3y[-1],l4y[-1])])
    #  # xch4 around box
    #  sq1=c(dat2[(ycent-i),l1x[-1]],dat2[l2y[-1],(xcent+i)],dat2[ceiling(ycent+i/4),l3x[-1]],dat2[l4y[-1],(xcent-50)])
    #}
    
    
    square=cbind(box,sq1)
    colnames(square)=c("x","y","xch4")
  }
  if(wdir1>135 & wdir1<225){
    # create box of side length i
    l1x=(xcent-i):(xcent+i)
    l1y=rep(ceiling(ycent-i/4),length(l1x))
    l2y=ceiling(ycent-i/4):(ycent+i)
    l2x=rep((xcent+i),length(l2y))
    l3x=(xcent+i):(xcent-i)
    l3y=rep((ycent+i),length(l3x))
    l4y=(ycent+i):ceiling(ycent-i/4)
    l4x=rep((xcent-i),length(l2y))
    # lat/lon coordinages of box
    box=cbind(xx[c(l1x[-1],l2x[-1],l3x[-1],l4x[-1])],yy[c(l1y[-1],l2y[-1],l3y[-1],l4y[-1])])
    # xch4 around box
    sq1=c(dat2[ceiling(ycent-i/4),l1x[-1]],dat2[l2y[-1],(xcent+i)],dat2[(ycent+i),l3x[-1]],dat2[l4y[-1],(xcent-i)])
    square=cbind(box,sq1)
    colnames(square)=c("x","y","xch4")
  }
  
  # component of u-wind in direction of square edge; 0 for horizontal lines, 1 for vertical lines
  u_fac=c(rep(0,length(l1x[-1])),rep(1,length(l2y[-1])),rep(0,length(l3x[-1])),rep(-1,length(l4y[-1])))
  # component of v-wind in direction of square edge; 1 for horizontal lines, 0 for vertical lines
  v_fac=c(rep(1,length(l1x[-1])),rep(0,length(l2y[-1])),rep(-1,length(l3x[-1])),rep(0,length(l4y[-1])))
  
  return( list(
    square=square,
    u_fac=u_fac,
    v_fac=v_fac)
  )
}

##-------------END function

## ------------ met.square --------------
# find met data at location of flux integral square 
met.square = function(met,met.nn,square) {
  ave.pbl=F
  if(ave.pbl==T){
    # Calculate level_h, the average height of the HRRR levels at the square where flux is calculated       
    if(length(dim(met$level_h0))==3){
      level_h=rep(NA,dim(met$level_h0)[3])
      for(k in 1:20){ #1st 20 levels will cover the PBL
        level_h[k]=mean(as.vector(met$level_h0[,,k])[met.nn])
      }
    }
    if(length(dim(met$level_h0))==2){
      level_h=rep(NA,dim(met$level_h0)[2])
      for(k in 1:20){
        level_h[k]=mean(as.vector(met$level_h0[,k])[met.nn])
      }
    }
    
    pblh=mean(as.vector(met$pblh0)[met.nn]) #mean PBL height at location of square
    
    levels=level_h[1:(max(which(level_h<pblh))+1)] # HRRR levels within PBL
    level_lim=0
    # Calculate boundaries between levels
    for(h in 2:length(levels)){
      level_lim=c(level_lim,level_lim[h-1]+2*(levels[h-1]-level_lim[h-1]))
    }
    z_level=diff(level_lim) # length of each level (m)
    
    # Calculate average wind speed within PBL, weighted by mass of air in each HRRR level
    uwind1=matrix(NA,length(met.nn),max(which(level_h<pblh))) # wind speed at each point along square (for nearest HRRR grid cell)
    vwind1=uwind1
    wspd1=uwind1
    n_bin=uwind1 #moles of air in HRRR level
    for(h in 1:max(which(level_h<pblh))){ #loop through HRRR levels that are less than PBL height
      if(length(dim(met$pres))==2){
        n_bin[,h]=met$pres[,h][met.nn]*100*z_level[h]/(8.314*met$temp[,h][met.nn])
        uwind1[,h]=met$uwind[,h][met.nn]*met$pres[,h][met.nn]*100*z_level[h]/(8.314*met$temp[,h][met.nn])
        vwind1[,h]=met$vwind[,h][met.nn]*met$pres[,h][met.nn]*100*z_level[h]/(8.314*met$temp[,h][met.nn])
        wspd1[,h]=met$wspd[,h][met.nn]*met$pres[,h][met.nn]*100*z_level[h]/(8.314*met$temp[,h][met.nn])
      }
      if(length(dim(met$pres))==3){
        n_bin[,h]=as.vector(met$pres[,,h])[met.nn]*100*z_level[h]/(8.314*as.vector(met$temp[,,h])[met.nn])
        uwind1[,h]=as.vector(met$uwind[,,h])[met.nn]*as.vector(met$pres[,,h])[met.nn]*100*z_level[h]/(8.314*as.vector(met$temp[,,h])[met.nn])
        vwind1[,h]=as.vector(met$vwind[,,h])[met.nn]*as.vector(met$pres[,,h])[met.nn]*100*z_level[h]/(8.314*as.vector(met$temp[,,h])[met.nn])
        wspd1[,h]=as.vector(met$wspd[,,h])[met.nn]*as.vector(met$pres[,,h])[met.nn]*100*z_level[h]/(8.314*as.vector(met$temp[,,h])[met.nn])
      }
    }
    n_PBL=apply(n_bin,1,sum) # mol air/m2 from ground to PBLH
    wspd2=apply(wspd1,1,sum)/n_PBL # ave wind speed from ground to PBLH at each point along square 
    
    # ave uwind and vwind for all locations on square
    uwind2=mean(apply(uwind1,1,sum)/n_PBL,na.rm=T) # ave uwind from HRRR from ground to PBLH at locations on square
    vwind2=mean(apply(vwind1,1,sum)/n_PBL,na.rm=T) # ave vwind from HRRR from ground to PBLH at locations on square
  }
  
  # moles of air through entire column at location of square     
  column_air=mean(as.vector(met$pres_surf)[met.nn]*100*1000/(9.80665*28.965))  #mol air/m2
  
  if(is.null(dim(met$uwind))){
    uwind2=met$uwind[3]
    vwind2=met$vwind[3]
  }else if(is.na(dim(met$uwind)[3])){
    uwind2=mean(as.vector(met$uwind[,3])[met.nn])
    vwind2=mean(as.vector(met$vwind[,3])[met.nn])
  }else{
    uwind2=mean(as.vector(met$uwind[,,3])[met.nn])
    vwind2=mean(as.vector(met$vwind[,,3])[met.nn])
  }
  
  surf_ht1=mean(as.vector(met$surf_ht)[met.nn])
  
  return( list(
    uwind2=uwind2,
    vwind2=vwind2,
    column_air=column_air,
    surf_ht=surf_ht1)
  )
}
## --------- END function --------------

## ------------ flux.calc --------------
# calculate flux divergence for individual square 
flux.calc = function(square.info, met.sq){
  # Calculate flux divergence around square; Eq. 6 from Conley et al. (2017) Atmos. Meas. Tech., https://doi.org/10.5194/amt-10-3345-2017
    # flux for single square of edge dx centered on cntr.x, cntr.y
    flux=(square.info$square[,"xch4"]-mean(square.info$square[,"xch4"],na.rm=T))*met.sq$column_air*met.sq$uwind2*square.info$u_fac*square.info$d1 + (square.info$square[,"xch4"]-mean(square.info$square[,"xch4"],na.rm=T))*met.sq$column_air*met.sq$vwind2*square.info$v_fac*square.info$d1
    flux0=sum(flux,na.rm=T)*16.04*60*60/1000 # convert mol/s to kg/hr
  return( flux0)
}
## --------- END function --------------



