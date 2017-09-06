# read in and plot ASASSN data
# alittle housekeeping
options(digits=12)
library("earth") # for MARS
library("crayon") # to add color to text

source("plot_funcs.R")
source("data_funcs.R")

# the name of the file
asassn.csv.file <- "data/ASASSN_latest.csv"
# cull limits to reject really wild points
mag.cull.limit = 0.05
# mag margin for expanding plot limits (positive number)
mag.margin <- 0.025
tmargin <- 1 # positive number, days extra on xscale. 
####### MARS
marsOrder <- 11 # maximum number of knots
marsPenalty <- 5 # set to 0 to avoid penalizing knots in pruning pass
#marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "backward" # set to "none" to avoid pruning
##########
earliestJD <- 2457449

source("input_files/VlineParams.R")
source("input_files/dip_mask.R")


#####################################################################
# load the files
asassn_data <- read.csv(asassn.csv.file,header=TRUE)
# clean obvious wild points
mag.median <- median(asassn_data$mag)
good_data <- (asassn_data$mag < (mag.median + abs(mag.cull.limit))) & (asassn_data$mag > (mag.median- abs(mag.cull.limit)))

#fit a spline after filtering out the dips
not.dips <- filter.dips.JD(asassn_data$HJD,dip.mask)
# set a nice round tmin
tmin <- max(min(asassn_data$HJD[good_data]),earliestJD)
tmin <- 20*floor(tmin/20)

desmat <- asassn_data$HJD[good_data & not.dips] - tmin
mars <- earth(x=desmat,y= asassn_data$mag[good_data & not.dips],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)

# plot limits
my.y.lims <- c(max(asassn_data$mag[good_data]) + mag.margin,min(asassn_data$mag[good_data]) - mag.margin)
my.xlab <- paste("Julian Date - ",tmin)
my.x.lims <- c(0,max(desmat) + tmargin)
# plot
quartz("asas-sn")
plot.times <- asassn_data$HJD[good_data] - tmin
plot(plot.times,asassn_data$mag[good_data],col="darkgreen",pch=20,main="ASASSN V Data",xlab=my.xlab,ylab="V magnitude",ylim=my.y.lims,xlim= my.x.lims)
# add grid
grid(col="black")

# plot MARS fir as a line
lines(x=desmat,y=mars$fitted.values,col= "black",lwd=2)
# plot vertical lines if any
jdLine <- jdLine - tmin
verticalDateLines(jdLine, jdLineText, my.y.lims, jdLineColor)
