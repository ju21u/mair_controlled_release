# MAIR Stanford Release Comparision 2022-02-21
#setwd(FindThisScriptsLocation())
cat('\014')
setwd("/Volumes/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/From_Steve/money_plot/")
#graphics.off()
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
	lm.york.all = YorkFit( X= Stanford[ , "cr_kgh_CH4_lower90" ] , 
	                       Y =  Stanford[ , "Emissions" ] , 
	                       Xstd = CI.cr / 2.3 ,
	                       Ystd = CI.mair  / 2.3 )
	
	
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
	Emiss.Ueff.Stanford = IME*S1[,"ueff.stanford"]
	Emiss.U.Stanford = IME*S1[,"U.stanford"]
	
	
	lm.sma = lmodel2( Emissions[lok] ~ cr_kgh_CH4_mean90[lok], data = Stanford )
	sma.coeffs= as.numeric( lm.sma$regression.results[3, 2:3] )
	
	lm.sma.all = lmodel2( Emissions ~ cr_kgh_CH4_mean90, data = Stanford )
	sma.coeffs.all = as.numeric( lm.sma.all$regression.results[3, 2:3] )
	


	## ===============  main figure ===============================
#pdf("Fig.Stanford_Mair_20220322.pdf", height=8.5, width=9.5,pointsize=18)

Hmisc::errbar(x= Stanford[, "cr_kgh_CH4_mean90"], y=Stanford[,"Emissions"],
	      yplus=Stanford[,"Upper"],yminus=Stanford[,"Lower"],pch=16,cex=1.3
              , xlab="Stanford Controlled Release (90 s mean)",ylab="MAIR IME Estimates ",
		 xlim=c(0,1000), ylim=c(0,1000), errbar.col = "darkgrey",font.lab=2,font.axis=2)
points( Stanford[ Stanford[, "Flag"] == "RM", c("cr_kgh_CH4_mean90","Emissions") ], col = "azure4", pch=15, cex=1.5)
points( Stanford[ l.nodetect , c("cr_kgh_CH4_mean90","Emissions") ], col = "red", lwd=3, pch = 15, cex = 1.5)
lines( c( t( horr.bars ) ), x3, col="darkgrey")
points(  c( t( horr.bars ) ), x3,col="darkgrey", pch="|")

title( main ="MAIR Controlled Release Validation")
#mtext(side=3,adj=.95,paste("MAIR Conf Int +/- =",round(Mair.err)))
mtext(side=3,adj=.03,expression(paste("Units: kg ",CH[4]," / hr")))
#points(Stanford[l.nodetect, c("cr_kgh_CH4_mean90", "Emissions")], pch=16,cex=.6,col="red")
#points(Stanford[l.both, c("cr_kgh_CH4_mean90", "Emissions")], pch=16,cex=.6,col="deepskyblue3")
## add symbols for those cases that include Gauss Integral
abline( 0, 1, col = "black", lty = 2, lwd = 3)

#abline( lm.york.all[1], lm.york.all[2], lty = 2, lwd = 3, col = "blue")
#abline( lm.york[1], lm.york[2], lty = 1, lwd = 3, col = "blue")
#abline( lm.york.2[1], lm.york.2[2], lty = 3)
abline( sma.coeffs.all[1], sma.coeffs.all[2], col="deepskyblue3",lwd = 3, lty = 2)
abline( sma.coeffs[1], sma.coeffs[2], col="deepskyblue3",lwd = 3)
legend( "bottomright", 
         legend = c("MAIR","Removed","No Detection & Removed","1:1",  "SMA Type II (filtered)", "SMA Type II (original)") ,
         pch=c(16, 15, 15, -1, -1,-1),
	 lty=c( -1, -1, -1, 2,  1, 2) ,
	 lwd = c(3, 3, 3, 3, 3, 3, 3),
	 #	 col=c("black","azure4","red", "black", "blue","deepskyblue3", "blue","deepskyblue3"))
	 col=c("black","azure4","red", "black","deepskyblue3","deepskyblue3"))
# 
# pch=c(16, 15, 15, -1, -1, -1,-1),
# lty=c( -1, -1, -1, 1, 3, 2,1) ,
# col=c("black","grey","red","blue","black","black","deepskyblue3"))
#text(800,850,"Slope\n=1")
#text(900, lm.york.2[1] + lm.york.2[2]*900,round(lm.york.2[2],2))
#text(900, lm.york[1] + lm.york[2]*950, round(lm.york[2],2), col = "blue")
text(890, 990, round(sma.coeffs[2],2),col="deepskyblue3")

#text(980, lm.york[1] + lm.york[2]*880, round(lm.york.all[2],2), col = "blue")
text(970, 920, round(sma.coeffs.all[2],2), col="deepskyblue3")


#dev.off()
#               ## end of main plot =====================================
# #  wind comnparison plot
# #pdf("Fig.winds_comp.pdf",width=13.5,height=8.5, pointsize=18)
# par(mfrow=c(1,2) )
# ##                   wind speed comparison
# plot(Stanford[,c("S.winds","LES.winds")],pch=16, cex=1.3,font.lab=2,font.axis=2)
#  abline(0,1)
#  text(3,3.2, "1:1" )
#  ##
#  ##                 Emissions difference vs wind difference
#  plot(X ,Y,
# 	xlab= expression(paste(Delta,"Wind (m/s, MAIR - Stanford)")),
# 	ylab= expression(paste(Delta,"Emissions ( kg/hr, MAIR - Stanford")),
# 				pch=16, cex=1.3,font.lab=2,font.axis=2)
# 	mtext( side=3,line=1, expression(paste(Delta,"Wind explains > 50% ",Delta,"Emission, r = .63" )) )
#  abline(my.coeff[1], my.coeff[2])
#  text(1.6,300, "SMA fit")
#  #dev.off()
# 
# 
all.flare = read.csv("/Volumes/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/Outputs/CSVs/20220326_summary.flare.1x1.sd.1.5.csv")
all.release = read.csv("/Volumes/GoogleDrive/My Drive/Research/Harvard_Research/MethaneAIR/Outputs/CSVs/20220326_summary.release.1x1.sd.1.5.csv")

names(all.flare) = names(all.release)
all.plumes = rbind(all.flare, all.release)
select_plumes = (all.plumes$Gaussian > 500 ) & (all.plumes$IME > 500 )  

lm.sma.DI_IME = lmodel2( Gaussian ~ IME, data = (all.plumes[select_plumes, ]) )
sma.coeffs= as.numeric( lm.sma.DI_IME$regression.results[3, 2:3] )

lm.sma.DI_IME = lmodel2( Gaussian ~ IME, data = (all.plumes) )
sma.coeffs= as.numeric( lm.sma.DI_IME$regression.results[3, 2:3] )



## ===============  main figure ===============================
#pdf("Fig.Stanford_Mair_20220322.pdf", height=8.5, width=9.5,pointsize=18)
all.plumes$Gaussian.u[is.na(all.plumes$Gaussian.l)] = NA
all.plumes$Gaussian.l[is.na(all.plumes$Gaussian.l)] = NA

horr.bars = cbind(all.plumes[, c( "IME.l","IME.u")], rep(NA, nrow(all.plumes)))
x3 = rep(all.plumes$Gaussian, each=3)  # the points around the horizontal bars will be plotted

Hmisc::errbar(x= all.plumes$IME, y=all.plumes$Gaussian,
              yplus= all.plumes$Gaussian.u, yminus=all.plumes$Gaussian.l, pch=16,cex=1.3
              , xlab="mIME",ylab="DI",
              xlim=c(0,2300), ylim=c(0,2300), errbar.col = "darkgrey",font.lab=2,font.axis=2)
points(all.release$IME, all.release$Gaussian, col = "azure4", pch=15, cex=1.5)
points( Stanford[ l.nodetect , c("cr_kgh_CH4_mean90","Emissions") ], col = "red", lwd=3, pch = 15, cex = 1.5)
lines( c( t( horr.bars ) ), x3, col="darkgrey")
points(  c( t( horr.bars ) ), x3,col="darkgrey", pch="|")

title( main ="DI vs mIME")
#mtext(side=3,adj=.95,paste("MAIR Conf Int +/- =",round(Mair.err)))
mtext(side=3,adj=.03,expression(paste("Units: kg ",CH[4]," / hr")))
#points(Stanford[l.nodetect, c("cr_kgh_CH4_mean90", "Emissions")], pch=16,cex=.6,col="red")
#points(Stanford[l.both, c("cr_kgh_CH4_mean90", "Emissions")], pch=16,cex=.6,col="deepskyblue3")
## add symbols for those cases that include Gauss Integral
abline( 0, 1, col = "black", lty = 2, lwd = 3)

#abline( lm.york.all[1], lm.york.all[2], lty = 2, lwd = 3, col = "blue")
#abline( lm.york[1], lm.york[2], lty = 1, lwd = 3, col = "blue")
#abline( lm.york.2[1], lm.york.2[2], lty = 3)
#abline( sma.coeffs.all[1], sma.coeffs.all[2], col="deepskyblue3",lwd = 3, lty = 2)
abline( sma.coeffs[1], sma.coeffs[2], col="deepskyblue3",lwd = 3)
legend( "bottomright", 
        legend = c("Flare","Release","No Detection & Removed","1:1", "SMA Type II (original)") ,
        pch=c(16, 15, 15, -1, -1),
        lty=c( -1, -1, -1, 2,  1) ,
        lwd = c(3, 3, 3, 3, 3, 3),
        #	 col=c("black","azure4","red", "black", "blue","deepskyblue3", "blue","deepskyblue3"))
        col=c("black","azure4","red", "black","deepskyblue3"))
# 
# pch=c(16, 15, 15, -1, -1, -1,-1),
# lty=c( -1, -1, -1, 1, 3, 2,1) ,
# col=c("black","grey","red","blue","black","black","deepskyblue3"))
#text(800,850,"Slope\n=1")
#text(900, lm.york.2[1] + lm.york.2[2]*900,round(lm.york.2[2],2))
#text(900, lm.york[1] + lm.york[2]*950, round(lm.york[2],2), col = "blue")
text(890, 990, round(sma.coeffs[2],2),col="deepskyblue3")

#text(980, lm.york[1] + lm.york[2]*880, round(lm.york.all[2],2), col = "blue")
#text(970, 920, round(sma.coeffs.all[2],2), col="deepskyblue3")


