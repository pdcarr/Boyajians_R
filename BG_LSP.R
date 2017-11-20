library(lomb)
library("MASS") # for rlm() and lqs()
source("special_funcs.R")

# lightcurve has been read in by the aavso plottng script with names autoassigned by read.csv()
# cleanBand is output fron the plotting script's cleaning process
fit.resids <- FALSE	# TRUE if doing the periodogram to the residuals vs. a fit
limit.period <- TRUE
max.period <- 30
plot.max.period <- 20 # days
set.ofac <- 8
use.frequency <- FALSE # set this to FALSE to plot vs. period
######## fitting parameters
bin.it <- FALSE
bg.n.knots <- 5

title.string <- "Bruce Gary Data L-S Periodogram"

# confined Gaussian window for binned data.
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
	our.resids <- (allSuperObs$V.mag - fit.points)
} else {
	desmat <- bg.data$MJD[ok.airmass] - tmin
	X <- bg.data$V.mag[ok.airmass]
#	theFit <- rlm(bg.data$V.mag[ok.airmass] ~ desmat,na.action="na.omit",psi=psi.bisquare)
    dipless <- filter.dips.JD(X,dip.mask)
    binWeights <- as.numeric(dipless) # weight of 1 if not in a known dip or flare, 0 otherwise.
	theFit <- smooth.spline(x=desmat,
						y= X,
						w= binWeights,
						all.knots=FALSE,
						nknots= bg.n.knots,
						keep.data=TRUE,
						cv=TRUE,
						penalty= 1)

	mytimes <- desmat
	fit.points <- predict(theFit,x= mytimes)$y
	our.resids <- (X - fit.points)

	# set up window for LSP
	N <- length(mytimes)
	w <- hamming.window(N)	# rectangular window
}
	
if (fit.resids) {

	quartz("residuals plot")
	plot(x= mytimes[dipless],
		y=our.resids[dipless],
		col='red', 
		pch=20,
		xlab="time",
		ylab="residual")
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
	mylsp <- lsp(x=X,
				times=mytimes,
				type="period",
				ofac=set.ofac,
				alpha=0.01,
				to=max.period)
	if(limit.period) {
		myxlims <- c(0, plot.max.period)
		plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string,xlim=myxlims)
	} else {
		plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string)
	}
}

grid(col="black")
