library(lomb)
# lightcurve has been read in by the aavso plottng script with names autoassigned by read.csv()
# cleanBand is output fron the plotting scrpts cleaning process: airmass, comparison stars, etc.
fit.resids <- TRUE
limit.period <- TRUE
max.period <- 40 # days
set.ofac <- 8
use.frequency <- FALSE # set this to FALSE to plot vs. period
desmat <- bg.data$MJD[ok.airmass] - tmin
theFit <- rlm(bg.data$V.mag[ok.airmass] ~ desmat,na.action="na.omit",psi=psi.bisquare)

mytimes <- desmat
quartz("LSP of Bruce Gary data")
title.string <- "Bruce Gary Data L-S Periodogram"

# confined Gaussian window
w <- cg.window(mytimes)
	
if (fit.resids) {
	X <- as.vector(theFit$residuals*t(w))
} else {
	X <- as.vector(bg.data$V.mag[ok.airmass]*t(w))
}


if(use.frequency) {
	mylsp <- lsp(x=theFit$residuals,times=mytimes,from=0.1,to=3.5,type="frequency",ofac=set.ofac,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Frequency (cycles/day)",ylab="normalized power",col="darkgreen",type="l",main= title.string)
} else {
	mylsp <- lsp(x=theFit$residuals,times=mytimes,type="period",ofac=set.ofac,alpha=0.01)
	if(limit.period) {
		myxlims <- c(0,max.period)
		plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string,xlim=myxlims)
	} else {
		plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string)
	}
}

grid(col="black")
