#### housekeeping
library("MASS") # for rlm() and lqs()
options(digits=12) # hard to read JDs without this setting
rm("superObs")
rm("allSuperObs")
source("input_files/VlineParams.R")
source("plot_funcs.R")

######################## control parameters
maxAirmass <- 2.0 # data with airmass higher than this will not be included
bin.width = 1/24 # days
t.epsilon = 1.0 # days
mag.epsilon <- 0.01 # magnitudes
bin.it <- TRUE
plot.col = "darkgreen"
plot.pch <- 20
dfile.name <- "data/BruceGary.csv"
############################################################

###### read in the data
bg.data <- read.csv(dfile.name,header=TRUE)

# filter out the high airmass data
ok.airmass <- bg.data$air.mass <= maxAirmass

# find earliest and latest times and highest and lowest magnitudes

# get the bin boundaries to be somewhere sensible

tmin <-  min(bg.data$MJD[ok.airmass],na.rm=TRUE)
if (floor(tmin*2) %% 2 == 0) {
    tmin <- floor(tmin*2 -1)/2
} else {
    tmin <- floor(tmin*2)/2
}

tmax = max(bg.data$MJD[ok.airmass],na.rm=TRUE)
x.limits <- c(0,(tmax-tmin) + t.epsilon)
y.limits <- c(max(bg.data$V.mag[ok.airmass],na.rm=TRUE) + mag.epsilon,
			 min(bg.data$V.mag[ok.airmass],na.rm=TRUE) - mag.epsilon)
x.label <- paste("MJD - ",tmin)

# if binning, do this:
if (bin.it) {
	plot.title <- paste("Bruce Gary Data binned with width",bin.width,"Days")
	y.label <- "Mean V over bin"
	allSuperObs <- data.frame(MJD=numeric(),V.mag=numeric(),Uncertainty=numeric(),stringsAsFactors=FALSE)
	superObs <- data.frame(MJD=numeric(),V.mag=numeric(),Uncertainty=numeric(),stringsAsFactors=FALSE)
	for (t in seq(tmin,tmax,by=bin.width)) {
		in.bin <- !is.na(bg.data$MJD) & bg.data$MJD >= t & bg.data$MJD < (t + bin.width) & ok.airmass
		n.in.bin <- length(bg.data$MJD[in.bin])
		if (n.in.bin > 0) {
			superObs[1,"MJD"] <- mean(bg.data$MJD[in.bin],na.rm=TRUE)
			superObs[1,"V.mag"] <- mean(bg.data$V.mag[in.bin],na.rm=TRUE)
			superObs[1,"Uncertainty"] <- mean(bg.data$SE[in.bin],na.rm=TRUE)/sqrt(n.in.bin)
			allSuperObs <- rbind(allSuperObs,superObs)
		} else {
			next
		}
	}
	# calculate robust linear fit
	desmat <- allSuperObs$MJD - tmin
	theFit <- rlm(allSuperObs$V.mag ~ desmat,na.action="na.omit",psi=psi.bisquare)
	quartz("binned Data Plot")
	plot(desmat,allSuperObs$V.mag,xlim=x.limits,ylim=y.limits,xlab=x.label,ylab=y.label,main=plot.title,pch= plot.pch,col=plot.col)
} else {
	quartz("Data Plot")
	y.label = "V Magnitude"
	plot.title <- paste("Bruce Gary Data with Airmass <=",maxAirmass)
	desmat <- bg.data$MJD[ok.airmass] - tmin
	theFit <- rlm(bg.data$V.mag[ok.airmass] ~ desmat,na.action="na.omit",psi=psi.bisquare)
	plot(desmat,bg.data$V.mag[ok.airmass],xlim=x.limits,ylim=y.limits,xlab=x.label,ylab=y.label,main=plot.title,pch= plot.pch,col=plot.col)
}

#### plot the fit 

myslope <- coefficients(theFit)[2]
basemag <- coefficients(theFit)[1]
curve(basemag + myslope*x,from=0,to=tmax-tmin, add=TRUE,col="black")
##### draw vertical lines at various dates
jdLine <- jdLine - 2400000.5 - tmin
if(drawDateLine) { verticalDateLines(jdLine, jdLineText, y.limits, jdLineColor)}

# add a grid
grid(col="black")

summary(theFit)
