########################### Parameters
#name of file
NW.file <- "data/NEOWISE_2017DR.csv"
# bin width in days
bin.width = 2 # days
plot.w1 = FALSE
w1.col = "blue"
w2.col="darkgreen"
bar.col="lightgrey"
pretty.mag <- 100
y.marg <- 0.02
########################### Process
#read the .csv file
neowise.data <- read.csv(file=NW.file,stringsAsFactors=FALSE)
#clean the data
ok.cc <- regexpr(pattern="^00",text= neowise.data$cc_flags,perl=TRUE) ==  1
# bin the w1 and w2 data

tmin = floor(min(neowise.data$mjd,na.rm=TRUE)) #floor time (modified Julian date)
tmax = ceil(max(neowise.data$mjd,na.rm=TRUE))	#ceiling time

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
y.limits <- c(ceil(max(allSuperObs$w2,na.rm=TRUE)*pretty.mag)/pretty.mag + y.marg,
				floor(min(allSuperObs$w2,na.rm=TRUE)*pretty.mag)/pretty.mag - y.marg)
t.limits <- c(tmin,tmax)
# plot the bins

w2.minus = allSuperObs$w2 - allSuperObs$w2.sigma
w2.plus = allSuperObs$w2 + allSuperObs$w2.sigma
quartz("NEOWISE binned")
plot.title <- "NEOWISE Data binned"

errbar(x=allSuperObs$mjd,
		y=allSuperObs$w2, 
		yplus=w2.plus,yminus=w2.minus,
		ylim=y.limits,xlim=t.limits,
		xlab= "MJD",ylab="W2 magnitude",
		main=plot.title,
		col=w2.col,errbar.col=bar.col,pch=20)
		
		
grid(col="black")