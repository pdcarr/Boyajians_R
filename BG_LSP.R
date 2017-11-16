library(lomb)
library("MASS") # for rlm() and lqs()
source("special_funcs.R")

# lightcurve has been read in by the aavso plottng script with names autoassigned by read.csv()
# cleanBand is output fron the plotting script's cleaning process
fit.resids <- TRUE
limit.period <- TRUE
max.period <- 20 # days
set.ofac <- 8
use.frequency <- FALSE # set this to FALSE to plot vs. period
######## fitting parameters
bg.n.knots <- 5

quartz("LSP of Bruce Gary data")
title.string <- "Bruce Gary Data L-S Periodogram"

# confined Gaussian window
if (bin.it) {
	desmat <- allSuperObs$MJD - tmin
	bg.binned.JDs <- allSuperObs$MJD  + 2400000.5
    dipless <- filter.dips.JD(bg.binned.JDs,dip.mask)
    binWeights <- as.numeric(dipless) # weight of 1 if not in a known dip or flare, 0 otherwise.
	theFit <- smooth.spline(x=desmat,
					y= allSuperObs$V.mag,
					w= binWeights,
					all.knots=FALSE,nknots= bg.n.knots,
					keep.data=TRUE,
					cv=TRUE,
					penalty= 1)
					
	cat("\n\n Smooth Spline Fit: \n")
	print(theFit $call)
	cat("\n")
	print(theFit $fit)
	fit.points <- predict(theFit,x=desmat)$y
	mytimes <- desmat
	w <- cg.window(mytimes)
	X <- allSuperObs$V.mag
} else {
	desmat <- bg.data$MJD[ok.airmass] - tmin
	theFit <- rlm(bg.data$V.mag[ok.airmass] ~ desmat,na.action="na.omit",psi=psi.bisquare)
	mytimes <- desmat
	w <- vector(n= length(mytimes),mode="numeric")
	w <- 1	# rectangular window
}
	
if (fit.resids) {
	X <- (allSuperObs$V.mag - fit.points)
	quartz("residuals plot")
	plot(x= mytimes[dipless],y=X[dipless],col='red', 
		pch=20,
		xlab="time",ylab="residual")
	points(x= mytimes[!dipless],y=X[!dipless],col='grey',pch=20)
	lines(x=mytimes,y=fit.points,col="darkgrey",lwd=2)
	grid(col="black")
	X <- as.vector(X*t(w))
} else {
	X <- as.vector(X*t(w))
}

quartz("Periodogram")
if(use.frequency) {
	mylsp <- lsp(x=X,times=mytimes,from=0.1,to=3.5,type="frequency",ofac=set.ofac,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Frequency (cycles/day)",ylab="normalized power",col="darkgreen",type="l",main= title.string)
} else {
	mylsp <- lsp(x=X,times=mytimes,type="period",ofac=set.ofac,alpha=0.01)
	if(limit.period) {
		myxlims <- c(0,max.period)
		plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string,xlim=myxlims)
	} else {
		plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string)
	}
}

grid(col="black")
