library("MASS") # for rlm() and lqs()
library("Hmisc")

########################### Parameters
#name of file
NW.file <- "data/NEOWISE_2017DR.csv"
# bin width in days
bin.width = 1 # days
plot.w1 = FALSE
w1.col = "blue"
w2.col="darkgreen"
bar.col="lightgrey"
pretty.mag <- 100
y.marg <- 0.02
worst.qual <- 10
########################### Process
#read the .csv file
neowise.data <- read.csv(file=NW.file,stringsAsFactors=FALSE)

#clean the data for contamination and confusion flags and low quality. flagged points won't be included in bins
ok.cc <- regexpr(pattern="^00",text= neowise.data$cc_flags,perl=TRUE) ==  1
ok.cc <- ok.cc & (neowise.data$qual_frame >= worst.qual)

# bin the w1 and w2 data

tmin = floor(min(neowise.data$mjd,na.rm=TRUE)) #floor time (modified Julian date)
tmax = ceiling(max(neowise.data$mjd,na.rm=TRUE))	#ceiling time

# create data frames
allSuperObs <- data.frame(mjd=numeric(),w1=numeric(),w1.sigma=numeric(),w2=numeric(),w2.sigma=numeric(),stringsAsFactors=FALSE)
superObs <- data.frame(mjd=numeric(),w1=numeric(),w1.sigma=numeric(),w2=numeric(),w2.sigma=numeric(),stringsAsFactors=FALSE)

# loop over all the times and average for each bin
for (t in seq(tmin,tmax,by=bin.width)) {
	in.bin <- !is.na(neowise.data$mjd) & neowise.data$mjd >= t & neowise.data$mjd < (t + bin.width) & ok.cc
	n.in.bin <- length(neowise.data$mjd[in.bin])
	if (n.in.bin > 0) {
		superObs[1,"mjd"] <- mean(neowise.data$mjd[in.bin],na.rm=TRUE)
		superObs[1,"w1"] <- mean(neowise.data$w1mpro[in.bin],na.rm=TRUE)
		superObs[1,"w1.sigma"] <- mean(neowise.data$w1sigmpro[in.bin],na.rm=TRUE)/sqrt(n.in.bin)
		superObs[1,"w2"] <- mean(neowise.data$w2mpro[in.bin],na.rm=TRUE)
		superObs[1,"w2.sigma"] <- mean(neowise.data$w2sigmpro[in.bin],na.rm=TRUE)/sqrt(n.in.bin)
		allSuperObs <- rbind(allSuperObs,superObs)
	} else {
		next
	}
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

errbar(x=allSuperObs$mjd,
		y=allSuperObs$w2, 
		yplus=w2.plus,yminus=w2.minus,
		ylim=y.limits,xlim=t.limits,
		xlab= "MJD",ylab="magnitude",
		main=plot.title,
		col=w2.col,errbar.col=bar.col,pch=20)
		
# put in a dashed line where the ALLWISE value is		
lines(x=neowise.data$mjd,y=neowise.data$w2mpro_allwise,lty="dashed",col=w2.col)

if(plot.w1) {
	points(x=allSuperObs$mjd,
			y=allSuperObs$w1,col=w1.col,pch=20)
	lines(x=neowise.data$mjd,y=neowise.data$w1mpro_allwise,lty="dotted",col=w2.col)
}

# add title string to top
title(plot.title)
grid(col="black")

quartz("W1 - W2")
plot(x=allSuperObs$mjd,
		y=allSuperObs$w1 - allSuperObs$w2, 
		xlim=t.limits,
		xlab= "MJD",ylab="magnitude difference",
		main="W1 - W2",
		col="darkred",pch=20)
		

grid(col="black")