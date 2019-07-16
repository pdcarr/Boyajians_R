library("MASS") # for rlm() and lqs()
library("Hmisc")
rm("plot.start.mjd")
rm("plot.stop.mjd")
########################### Parameters
#name of file
NW.file <- "data/NEOWISE_2019DR.csv"
# bin width in days
bin.width <- 24 # days
plot.w1 = TRUE
w1.col = "blue"
w2.col="green"
bar.col="lightgrey"
w1.pch <- 0 #symbold for plotting W1
w2.pch <- 1 # symbol for plotting W2
w2.offset <- 0.1
chi2.thresh <- 4 # only plot points with chi2 lower than this threshold
pretty.mag <- 100
y.marg <- 0.05
worst.qual <- 10
plot.start.mjd <- NA
plot.stop.mjd <- NA
#plot.start.mjd <- 57884
#plot.stop.mjd <- 57894
plot.diff <- FALSE # TRUE if you want to plot W1 - W2
source("input_files/VlineParams.R")
########################### Process
#read the .csv file
neowise.data <- read.csv(file=NW.file,stringsAsFactors=FALSE,header = TRUE)

#clean the data for contamination and confusion flags and low quality. flagged points won't be included in bins
ok.cc <- regexpr(pattern="^00",text= neowise.data$cc_flags,perl=TRUE) ==  1
ok.cc <- ok.cc & (neowise.data$qual_frame >= worst.qual)


tmin = floor(min(neowise.data$mjd,na.rm=TRUE)) #floor time (modified Julian date)
if(!is.na(plot.start.mjd)) {tmin <- max(tmin,plot.start.mjd)}
tmax = ceiling(max(neowise.data$mjd,na.rm=TRUE))	#ceiling time
if(!is.na(plot.stop.mjd)) {tmax <- min(tmax,plot.stop.mjd)}
# bin the w1 and w2 data

# create data frames
allSuperObs <- data.frame(mjd=numeric(),w1=numeric(),w1.sigma=numeric(),w2=numeric(),w2.sigma=numeric(),nobs.w1=numeric(),nobs.w2=numeric(),stringsAsFactors=FALSE)
superObs <- data.frame(mjd=numeric(),w1=numeric(),w1.sigma=numeric(),w2=numeric(),w2.sigma=numeric(),nobs.w1=numeric(),nobs.w2=numeric(),stringsAsFactors=FALSE)

# loop over all the times and average for each bin for W2
good.chi2.2 <- neowise.data$w2rchi2 <= chi2.thresh # check the Chi^2 threhold for W2
good.chi2.1 <- neowise.data$w1rchi2 <= chi2.thresh # check the Chi^s threhold for W1
for (t in seq(tmin,tmax,by=bin.width)) {
	in.bin <- !is.na(neowise.data$mjd) & neowise.data$mjd >= t & neowise.data$mjd < (t + bin.width) & ok.cc
#	n.in.bin <- length(neowise.data$mjd[in.bin])
	if (sum(in.bin & (good.chi2.1 | good.chi2.2)) > 0) {
		superObs[1,"mjd"] <- mean(neowise.data$mjd[in.bin],na.rm=TRUE) # there is at least 1 good data point in this bin
	} else {
		next # no good data points in this bin
	}
	if(sum(in.bin & good.chi2.1) > 0) {
		superObs[1,"nobs.w1"] <- sum(in.bin & good.chi2.1)
		superObs[1,"w1"] <- mean(neowise.data$w1mpro[in.bin & good.chi2.1],na.rm=TRUE)
		superObs[1,"w1.sigma"] <- mean(neowise.data$w1sigmpro[in.bin & good.chi2.1],na.rm=TRUE)/sqrt(sum(in.bin & good.chi2.1))
	} else {
		superObs[1,"nobs.w1"] <- 0
		superObs[1,"w1"] <- NaN
		superObs[1,"w1.sigma"] <- NaN
	}
	if(sum(in.bin & good.chi2.2) > 0) {
		superObs[1,"nobs.w2"] <- sum(in.bin & good.chi2.2)
		superObs[1,"w2"] <- mean(neowise.data$w2mpro[in.bin & good.chi2.2],na.rm=TRUE) + w2.offset
		superObs[1,"w2.sigma"] <- mean(neowise.data$w2sigmpro[in.bin & good.chi2.2],na.rm=TRUE)/sqrt(sum(in.bin & good.chi2.2))
	} else {
		superObs[1,"nobs.w2"] <- 0
		superObs[1,"w2"] <- NaN
		superObs[1,"w2.sigma"] <- NaN
	}
	allSuperObs <- rbind(allSuperObs,superObs)

}

# calculate x and y limits
y.limits <- c(ceiling(max(c(allSuperObs$w1,allSuperObs$w2),na.rm=TRUE)*pretty.mag)/pretty.mag + y.marg,
				floor(min(c(allSuperObs$w1,allSuperObs$w2),na.rm=TRUE)*pretty.mag)/pretty.mag - y.marg)
t.limits <- c(tmin,tmax)

# plot the bins

w2.minus = allSuperObs$w2 - allSuperObs$w2.sigma
w2.plus = allSuperObs$w2 + allSuperObs$w2.sigma
quartz("NEOWISE binned")

if(plot.w1) {
	band.string <- "W1 & W2"
}else {
	band.string <- "W2"
}

plot.title <- paste("NEOWISE",band.string,"Data",bin.width,"day bins")
if (w2.offset != 0) {
	w2.leg.string <- paste("W2 + ",as.character(w2.offset),sep="")
} else {
	w2.leg.string <- "W2"
}
errbar(x=allSuperObs$mjd,
		y=allSuperObs$w2, 
		yplus=w2.plus,yminus=w2.minus,
		ylim=y.limits,xlim=t.limits,
		xlab= "MJD",ylab="magnitude",
		main=plot.title,
		col=w2.col,errbar.col=bar.col,pch=w2.pch)
		
# put in a dashed line where the ALLWISE value is		
lines(x=neowise.data$mjd,y=neowise.data$w2mpro_allwise +  w2.offset,type="l",lwd=1,lty="dashed",col=w2.col)

if(plot.w1) {
	w1.plus <- allSuperObs$w1 + allSuperObs$w1.sigma
	w1.minus <- allSuperObs$w1 - allSuperObs$w1.sigma	
	errbar(x=allSuperObs$mjd,
		y=allSuperObs$w1, 
		yplus=w1.plus,yminus=w1.minus,
		ylim=y.limits,xlim=t.limits,
		col=w1.col,errbar.col=bar.col,pch=w1.pch,add=TRUE)

	# points(x=allSuperObs$mjd,
			# y=allSuperObs$w1,col=w1.col,pch=20)
	lines(x=neowise.data$mjd,y=neowise.data$w1mpro_allwise,type="l",lwd=1,lty="dashed",col=w1.col)
}
# add vertical event lines
##### draw vertical lines at various dates
if(drawDateLine) { 
	my.vlines <- jdLine - 2400000.5
	verticalDateLines(my.vlines, jdLineText, y.limits, jdLineColor)
	}


# add title string to top
title(plot.title)
grid(col="black")

# legend

legend(x= t.limits[1], y = y.limits[2], legend=c("W1",w2.leg.string),col = c(w1.col,w2.col),
       border = "black", pch=c(w1.pch,w2.pch),cex = 1)

if(plot.diff) {
	quartz("W1 - W2")
	plot(x=allSuperObs$mjd,
			y=allSuperObs$w1 - allSuperObs$w2, 
			xlim=t.limits,
			xlab= "MJD",ylab="magnitude difference",
			main="W1 - W2",
			col="darkred",pch=20)
			
	
	grid(col="black")
}