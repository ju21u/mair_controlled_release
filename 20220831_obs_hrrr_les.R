# MAIR Stanford Release Comparision 2022-02-21
#setwd(FindThisScriptsLocation())
cat('\014')
setwd("/Volumes/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/money_plot/")
graphics.off()
library(lubridate)
library(Hmisc)
library(lmodel2)
source("York.R")
# read Stanford file
	stanford = read.table("matchedDF_MAIR_unblindedToMAIR_2022-03-23.csv",header=T,sep=",")
#  this file combines RF04_flagged.csv, RF05_flagged.csv; check that it is correct
	mair = read.table("MAIR_release_data_2022-03-23.csv",header=T,sep=",")
#

	#convert dates in character format to POSIXct
		t0 = as_datetime( stanford[,"Stanford_timestamp"] )
	# convert Ju's short-form dates
	my.day=rep( "2021-08-03 ", nrow(mair))
	my.day[ 1:9] = "2021-07-30 "
	my.day.time = paste(my.day, mair[,"Time"], ":00" , sep = "")
		tt = as_datetime( my.day.time )

	## combind the data frames ------------------------------------------------------
Stanford = data.frame( stanford, mair )
		# note-- need to check this step since time matching was not straightforward
		# they are being matched by segment
	## combind the data frames ------------------------------------------------------
#  instert the POSIXct times
	Stanford[,"Stanford_timestamp"] = t0
	Stanford[,"Time"] = tt
# confidence interval vector
	CI.cr = ( - Stanford[, "cr_kgh_CH4_lower90"] + Stanford[,"cr_kgh_CH4_upper90"])
	CI.mair = ( - Stanford[ , "Lower"] + Stanford[,"Upper"])
#  change no-detect from NA to 0 with a confidence interval
	       l.nodetect = is.na( Stanford[,"Emissions"])
        Mair.err = mean((Stanford[,"Upper"]-Stanford[,"Lower"])/2,na.rm=T)
        Stanford[ l.nodetect, c("Emissions", "Lower", "Upper")] = c(0, 0, Mair.err)
#  make a matrix of the Stanford CIs, add a column of NAs at the right;  for plotting
	horr.bars = cbind(Stanford[, c( "cr_kgh_CH4_lower90","cr_kgh_CH4_upper90")],rep(NA,nrow(Stanford)))
	x3 = rep( mair[,"Emissions"],each=3)  # the points around the horizontal bars will be plotted

	lok = ! Stanford[,"Flag"] == "RM" & ( ! l.nodetect )
	lm.york = YorkFit( X= Stanford[ lok, "cr_kgh_CH4_lower90" ] , 
			  Y =  Stanford[ lok, "Emissions" ] , 
			  Xstd = CI.cr[ lok ] / 2.3 ,
			  Ystd = CI.mair[ lok ] / 2.3 )

	lm.york.2 = YorkFit( X= Stanford[ , "cr_kgh_CH4_lower90" ] , 
			  Y =  Stanford[ , "Emissions" ] , 
			  Xstd = CI.cr[ ] / 2.3 ,
			  Ystd = CI.mair[ ] / 2.3 )

 	X = -1* (Stanford[,"S.winds"]-Stanford[,"LES.winds"] )
 	Y = -1* (Stanford[,"cr_kgh_CH4_mean90"] - Stanford[,"Emissions"] )
 	#X2 = -1* (Stanford[,"S.winds"]-Stanford[,"LES.winds"] )/Stanford[ , "S.winds"]
 	#Y2 = -1* (Stanford[,"cr_kgh_CH4_mean90"] - Stanford[,"Emissions"] )/ Stanford[, "cr_kgh_CH4_mean90"]
	lm.2 = lmodel2 ( Y ~ X )
	my.coeff = as.numeric( lm.2$regression.results[3, 2:3 ] )
	#lm.2.2 = lmodel2 ( Y2 ~ X2 )
	#my.coeff.2 = as.numeric( lm.2.2$regression.results[3, 2:3 ] )

	# read the winds and add to the data frame
	xx=read.table("summary.release.winds.csv",header=T,sep=",")
	S1=cbind(Stanford,xx)
	colnames(S1)[37]="U.stanford"
	# compute the IME since it was not in Ju's data frame
	IME=S1[,"Emissions"]/S1[,"ueff.LES"]
	
	Emiss.Ueff.HRRR = IME*S1[,"ueff.HRRR"]
	Emiss.U.HRRR = IME*S1[,"HRRR"]
	
	Emiss.U.LES = IME*S1[,"WS"]
	
	
	Emiss.Ueff.Stanford = IME*S1[,"ueff.stanford"]
	Emiss.U.Stanford = IME*S1[,"U.stanford"]
	lm.sma = lmodel2( Emissions ~ cr_kgh_CH4_mean90, data = Stanford )
	sma.coeffs= as.numeric( lm.sma$regression.results[2, 2:3] )



	## ===============  main figure ===============================
#pdf("Fig.Stanford_Mair_20220322.pdf", height=8.5, width=9.5,pointsize=18)
# 
# Hmisc::errbar(x= Stanford[, "cr_kgh_CH4_mean90"], y=Stanford[,"Emissions"],
# 	      yplus=Stanford[,"Upper"],yminus=Stanford[,"Lower"],pch=16,cex=1.3
#               , xlab="Stanford Controlled Release (90 s mean)",ylab="MAIR IME + Emission",
# 		 xlim=c(0,1000), ylim=c(0,1000), errbar.col = "darkgrey",font.lab=2,font.axis=2)
# points( Stanford[ Stanford[, "Flag"] == "RM", c("cr_kgh_CH4_mean90","Emissions") ], col = "grey",pch=15,cex=1.12)
# points( Stanford[ l.nodetect , c("cr_kgh_CH4_mean90","Emissions") ], col = "green",lwd=3)
# lines( c( t( horr.bars ) ), x3, col="darkgrey")
# points(  c( t( horr.bars ) ), x3,col="darkgrey", pch="|")
# 
# title( main ="MAIR Controlled Release Validation")
# #mtext(side=3,adj=.95,paste("MAIR Conf Int +/- =",round(Mair.err)))
# #mtext(side=3,adj=.05,expression(paste("Units: kg ",CH[4]," / hr")))
# points(Stanford[l.nodetect, c("cr_kgh_CH4_mean90", "Emissions")], pch=16,cex=.6,col="cyan")
# #points(Stanford[l.both, c("cr_kgh_CH4_mean90", "Emissions")], pch=16,cex=.6,col="orange")
# ## add symbols for those cases that include Gauss Integral
# abline( 0, 1 )
# abline( lm.york[1], lm.york[2], lty = 2)
# abline( lm.york.2[1], lm.york.2[2], lty = 3)
# abline( sma.coeffs[1], sma.coeffs[2], col="orange",lwd=2.5)
# legend( "bottomright", 
#          legend = c("MAIR","MAIR ROM","MAIR No Detect","1:1","York Fit, all", "York, Endpoint > 2","SMA Type II Regression") ,
#          pch=c(16, 15, 16, -1, -1, -1,-1), 
# 	 lty=c( -1, -1, -1, 1, 3, 2,1) ,
# 	 col=c("black","grey","green","black","black","black","orange"))
# text(800,850,"Slope\n=1")
# text(900, lm.york.2[1] + lm.york.2[2]*900,round(lm.york.2[2],2))
# text(950, lm.york[1] + lm.york[2]*950,round(lm.york[2],2))
# text(913,968,round(sma.coeffs[2],2),col="orange")
#dev.off()
              ## end of main plot =====================================
#  wind comnparison plot
#pdf("Fig.winds_comp.pdf",width=13.5,height=8.5, pointsize=18)
par(mfrow=c(1,2) )
##                   wind speed comparison
plot(Stanford[,c("S.winds","LES.winds")],pch=16, cex=1.3,
     font.lab=2,font.axis=2,
     xlab = expression("Stanford winds (m/s)"),
     ylab= expression("Input winds")
                      )
points(Stanford[,c("S.winds","HRRR.winds")], pch=18, cex=1.9, col = "red")
legend( "topright", 
        legend = c("WRF-LES-HRRR","HRRR") ,
                 pch=c(16, 18 ),
        	 lty=c( -1, -1) ,
        	 col=c("black","red"))
 abline(0,1)
 text(3,3.2, "1:1" )
 

 ##
 ##                 Emissions difference vs wind difference
 plot(X ,Y,
	xlab= expression(paste(Delta,"Wind (m/s, MAIR - Stanford)")),
	ylab= expression(paste(Delta,"Emissions ( kg/hr, MAIR - Stanford")),
				pch=16, cex=1.3,font.lab=2,font.axis=2)
	mtext( side=3,line=1, expression(paste(Delta,"Wind explains > 50% ",Delta,"Emission, r = .63" )) )
 abline(my.coeff[1], my.coeff[2], col = "black")
 text(1.6,300, "SMA fit")
 #dev.off()



 ######### EXPERIMENTAL  ######### 
 ##                 Emissions difference vs wind difference
 
 
 X = -1* (Stanford[,"S.winds"]-Stanford[,"LES.winds"] )
 Y = -1* (Stanford[,"cr_kgh_CH4_mean90"] - Emiss.Ueff.Stanford )
 lm.2 = lmodel2 ( Y ~ X )
 my.coeff = as.numeric( lm.2$regression.results[3, 2:3 ] )
 
 points(X ,Y,
      xlab= expression(paste(Delta,"Wind (m/s, MAIR - Stanford)")),
      ylab= expression(paste(Delta,"Emissions ( kg/hr, MAIR - Stanford")),
      pch=18, cex=1.8,font.lab=2,font.axis=2, col = "red")
 mtext( side=3,line=1, expression(paste(Delta,"Wind explains > 50% ",Delta,"Emission, r = .63" )) )
 abline(my.coeff[1], my.coeff[2], col = "red")
 #abline(0, 1, col = "grey", lty = 2)
 text(1.6,300, "SMA fit")
 
 
 
 ######### EXPERIMENTAL 2   ######### 
 ##  Stanford Emiss vs IME-HRRR Ueff
 ##  Stanford Emiss vs IME-WRF-LES-HRRR Ueff
 ##  Stanford Emiss vs Stanford wind Ueff
 
 ##  Stanford Emiss vs IME-HRRR U
 ##  Stanford Emiss vs IME-WRF-LES-HRRR U
 ##  Stanford Emiss vs Stanford wind U
 
 X = Stanford[,"cr_kgh_CH4_mean90"]
 
 Y = Stanford[,"Emissions"]
 Y2 = Emiss.Ueff.HRRR 
 Y3 = Emiss.Ueff.Stanford
 Y3[17] = NA
 Y4 = Emiss.U.LES
 Y5 = Emiss.U.HRRR
 Y6 = Emiss.U.Stanford
 Y6[17] = NA
 
 case_names =  c( "WRF-LES-HRRR Winds", "HRRR Winds", "Stanford Winds") #c("LES ueff", "HRRR ueff", "Stanford ueff", "LES u", "HRRR u", "Stanford u")
 
 par(mfrow=c(1, 3) )
 current.Y = Y
 lm.2 = lmodel2 ( current.Y ~ X )
 my.coeff = as.numeric( lm.2$regression.results[3, 2:3 ] )
 plot(X ,current.Y,
      xlab= "Stanford Reported Emission (kg/hr)",
      ylab= "Estimated Emissions (kg/hr)" ,
      xlim = c(0, 800),
      ylim = c(0, 800),
      pch=18, cex=1.8,font.lab=2,font.axis=2, col = "red")
 mtext( side=3,line=1, bquote( .(case_names[1])* " "*  ( R^2 == .(round(lm.2$rsquare, 3)) ) ) )
#  mtext( side=3,line=1, expression(paste( case_names[1], "(r2 =", round(lm.2$rsquare, 3), ")" )) )
 abline(my.coeff[1], my.coeff[2], col = "red")
 text(500, 300, "SMA fit")
 text(500,200,  bquote( "Slope" == .(round(lm.2$regression.results[3, "Slope"], 3) ) ))
 
 for (current.idx in seq(2, 3)) {
   current.Y = get(paste("Y", current.idx, sep = "")) 
   lm.2 = lmodel2 ( current.Y ~ X )
   my.coeff = as.numeric( lm.2$regression.results[3, 2:3 ] )
   plot(X ,current.Y,
        xlab= "Stanford Reported Emissions (kg/hr)",
        ylab= "Estimated Emissions (kg/hr)" ,
        xlim = c(0, 800),
        ylim = c(0, 800),
          pch=18, cex=1.8,font.lab=2,font.axis=2, col = "red")
  mtext( side=3,line=1, bquote( .(case_names[current.idx]) * " "* ( R^2 == .(round(lm.2$rsquare, 3)) ) ) )
   #mtext( side=3,line=1, (paste(case_names[current.idx], "($r^2$ =", round(lm.2$rsquare, 3), ")" )) )
   abline(my.coeff[1], my.coeff[2], col = "red")
   text(500,300, "SMA fit")
   text(500,200,  bquote( "Slope" == .(round(lm.2$regression.results[3, "Slope"], 3) ) ))
 }
 
 
 # TEST THE DIFFERENCE 
 
 # Stanford vs LES-EMISSION
 t.test(X, Y, paired = T)
 t.test(X, Y2, paired = T)
 t.test(X, Y3, paired = T)
 
 