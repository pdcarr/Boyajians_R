library(lomb)
# lightcurve has been read in by the aavso plottng script with names autoassigned by read.csv()
# cleanBand is output fron the plotting scrpts cleaning process: airmass, comparison stars, etc.

use.frequency <- TRUE # set this to FALSE to plot vs. period
mytimes <- lightcurve$JD[cleanBand[,1]] - tmin

title.string <- paste("AAVSO",allBands$bandinQ[1],"Data L-S Periodogram",sep=" ")
if(use.frequency) {
	mylsp <- lsp(x=thisFit$residuals,times=mytimes,from=0.1,to=3.5,type="frequency",ofac=4,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Frequency (cycles/day)",ylab="normalized power",col=allBands$plotColor[1],type="l",main= title.string)
} else {
	mylsp <- lsp(x=thisFit$residuals,times=mytimes,type="period",ofac=4,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Frequency (cycles/day)",ylab="normalized power",col=allBands$plotColor[1],type="l",main=title.string)

}
grid(col="black")
