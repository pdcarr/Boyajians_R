library(lomb)
# lightcurve has been read in by the aavso plottng script with names autoassigned by read.csv()
# cleanBand is output fron the plotting scrpts cleaning process: airmass, comparison stars, etc.

use.frequency <- FALSE # set this to FALSE to plot vs. period
desmat <- bg.data$MJD[ok.airmass] - tmin
theFit <- rlm(bg.data$V.mag[ok.airmass] ~ desmat,na.action="na.omit",psi=psi.bisquare)

mytimes <- desmat

title.string <- "Bruce Gary Data L-S Periodogram"
if(use.frequency) {
	mylsp <- lsp(x=theFit$residuals,times=mytimes,from=0.1,to=3.5,type="frequency",ofac=4,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Frequency (cycles/day)",ylab="normalized power",col=allBands$plotColor[1],type="l",main= title.string)
} else {
	mylsp <- lsp(x=theFit$residuals,times=mytimes,type="period",ofac=4,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="normalized power",col=allBands$plotColor[1],type="l",main=title.string)

}
grid(col="black")
