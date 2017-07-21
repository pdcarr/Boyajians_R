# this is designed to plot AAVSO magnitudes after some data cleaning, with decreasing magnitude, and to fit the data to a straight line.

# housekeeping
# options
options(digits=12) # hard to read JDs without this setting
# load the required packages
#library("MASS") # for rlm() and lqs()
library("earth") # for MARS
#library("crayon") # to add color to text


##################
# input parameters
mands.curve.name <- "data/Montet_SIMON_FFI_table.csv"

plotRelTimes <- TRUE
#####################
mands.JD.base <- 2454833
startPlot <- 0
########################

####### MARS
marsOrder <- 15
marsPenalty <- 4 # set to 0 to avoid penalizing knots in pruning pass
marsPMethod <- "none" # set to "none" to avoid pruning#marsPMethod <- "backward" # set to "none" to avoid pruning
marsPMethod <- "backward" # set to "none" to avoid pruning#marsPMethod <- "backward" # set to "none" to avoid pruning
##############################
### option to draw a vertical date line
drawDateLine <-  TRUE
jdLine <- c(1206,1520)
jdLineColor <- "red"
jdLineText <- c("D1200","D1520")
##########################################################################


# read in Montet and Simon FFI data:
MandSFFI <- read.csv(file= mands.curve.name,header=TRUE)



# create a vector of shifted times
#MandS.shift.times <- MandSFFI$Time + mands.JD.base
fit.times = MandSFFI$Time 
MandS.Fit <- earth(x= fit.times,y= MandSFFI$Normalized.Flux,nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)

print(summary(MandS.Fit))

##################################### Plot This Stuff ##################################
# axis limits

myxlims = c(startPlot,max(fit.times,na.rm=TRUE))

# calculate pretty y limits

myYlims = c(floor(min(MandSFFI$Normalized.Flux,na.rm=TRUE)*100)/100,ceiling(100*max(MandSFFI$Normalized.Flux,na.rm=TRUE))/100) # set up Y limits for reversed Y axis

# set up plot title text
myPlotTitle <- "Kepler FFI data for KIC 8462852 per Montet and Simon"


plot(fit.times, MandSFFI$Normalized.Flux,col="darkgreen",xlab=paste("Days from",mands.JD.base),ylab="Normalized Flux",xlim= myxlims,ylim = myYlims,main=myPlotTitle,pch=3,cex.main=0.7)
lines(fit.times,MandS.Fit$fitted.values,lwd=2,col="black")
grid(col="black")	

if(drawDateLine) {
	# draw a vertical line for a date of interest
	for(jd.index in 1:length(jdLine)) {
		lines(x=c(jdLine[jd.index],jdLine[jd.index]),y=myYlims,col=jdLineColor,lwd=1,lty="dashed")
		text(x= jdLine[jd.index],y=myYlims[2],labels=jdLineText[jd.index],pos=3,cex=0.5)
	}
}



