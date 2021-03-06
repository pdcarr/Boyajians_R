#### housekeeping
library("MASS") # for rlm() and lqs()
library("earth")
options(digits=12) # hard to read JDs without this setting
rm("superObs")
rm("allSuperObs")
source("input_files/VlineParams.R")
source("plot_funcs.R")

######################## control parameters
maxAirmass <- 2.0 # data with airmass higher than this will not be included
bin.width = 1 # days
#bin.width <- 15/1440 # days
t.epsilon = 1.0 # days
mag.epsilon <- 0.01 # magnitudes
bin.it <- TRUE
plot.col = "darkgreen"
plot.pch <- 20
#dfile.name <- "data/BruceGary.csv"
dfile.name <- "data/BG_gprime_latest.txt"
data.type <- "G prime"
#data.type <- "V"
source("input_files/dip_mask.R")
########## MARS control
n.knots <- 21
knot.penalty <- 6
min.span <- 14
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
y.limits <- c(max(bg.data$V.mag[ok.airmass],na.rm=TRUE) - mag.epsilon,
			 min(bg.data$V.mag[ok.airmass],na.rm=TRUE) + mag.epsilon)
x.label <- paste("MJD - ",tmin)

# if binning, do this:
if (bin.it) {
	plot.title <- paste("Bruce Gary Data binned with width",bin.width,"Days")
	y.label <- paste("Mean",data.type,"over bin")
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
	# mask dips from the fit
	bg.binned.JDs <- allSuperObs$MJD  + 2400000.5
    dipless <- filter.dips.JD(bg.binned.JDs,dip.mask)
    binWeights <- as.numeric(dipless) # weight of 1 if not in a known dip, 0 otherwise.
	# calculate robust linear fit
	desmat <- allSuperObs$MJD - tmin
    #theFit <- rlm(allSuperObs$V.mag ~ desmat,na.action="na.omit",psi=psi.bisquare,subset=dipless)
    theFit <- earth(x=desmat,
                    y=allSuperObs$V.mag,
                    weights=binWeights,
                    pmethod="exhaustive",
                    nk=n.knots,
                    penalty=knot.penalty,
                    minspan=min.span,
                    thresh=0.0001)

	# let's see what the fit is
	print(summary(theFit))
	# open a window to plot in
	quartz("binned Data Plot")
	# plot bins used in the fit
	plot(desmat[dipless],allSuperObs$V.mag[dipless],
		xlim=x.limits,ylim=y.limits,
		xlab=x.label,ylab=y.label,
		main=plot.title,
		pch= plot.pch,col=plot.col)
	# plot omitted bins
	points(x=desmat[!dipless],
		y=allSuperObs$V.mag[!dipless],
		pch= plot.pch,col="grey")
# not binned
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
#curve(basemag + myslope*x,from=0,to=tmax-tmin, add=TRUE,col="black")
lines(x=desmat,y=theFit$fitted.values,col="black")
##### draw vertical lines at various dates
jdLine <- jdLine - 2400000.5 - tmin
if(drawDateLine) { verticalDateLines(jdLine, jdLineText, y.limits, jdLineColor)}

# add a grid
grid(col="black")

summary(theFit)
