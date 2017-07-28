library(lomb)

# Kepler lightcurve (klightcurve) has been read in and cleaned up by the plottng script with names assigned
# cleanBand is output fron the plotting scrpts cleaning process: airmass, comparison stars, etc.

savdig <- options("digits")$digits
options(digits = 7)
use.frequency <- FALSE # set this to FALSEto plot vs. period

notthese <- is.na(klightcurve$Time)

# filter out data from any identitified dips

if (!is.null(dips)) {
    for (index in 1:length(dips$dstart)) {
        notthese <-  notthese & (klightcurve$Time >= dips$dstart[index] & klightcurve$Time <= dips$dstop[index])
    }

    in.bounds <- in.bounds & !notthese
}

# do the LS periodogram

mytimes <- klightcurve$Time[in.bounds]
t <- c(min(klightcurve$Time[in.bounds],na.rm=TRUE),max(klightcurve$Time[in.bounds],na.rm=TRUE))

title.string <- paste(c("Kepler Data L-S Periodogram",paste("From D",t[1]," to D",t[2],sep="")),collapse="\n")

if(use.frequency) {
	mylsp <- lsp(x=klightcurve$PDCSAP_flux[in.bounds],times=mytimes,from=0.1,to=3.5,type="frequency",ofac=4,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Frequency (cycles/day)",ylab="normalized power",col="darkgreen",type="l",main= title.string)
} else {
	mylsp <- lsp(x= klightcurve$PDCSAP_flux[in.bounds],times=mytimes,type="period",ofac=4,alpha=0.01)
	plot(mylsp$scanned,mylsp$power,xlab="Period (Days)",ylab="normalized power",col="darkblue",type="l",main=title.string)

}
grid(col="black")
options(digits=savdig)
