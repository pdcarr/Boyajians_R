# read in and plot ASASSN data
# alittle housekeeping
options(digits=12)
library("earth") # for MARS
library("crayon") # to add color to text

# the name of the file
asassn.csv.file <- "data/light_curve_ab426ded-039a-4c44-8976-ef7fa220899f.csv"
# cull limits to reject really wild points
mag.cull.limit = 0.05
# mag margin for expanding plot limits (positive number)
mag.margin <- 0.025
####### MARS
marsOrder <- 21 # maximum number of knots
marsPenalty <- 4 # set to 0 to avoid penalizing knots in pruning pass
marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "backward" # set to "none" to avoid pruning
#####################################################################
# load the files
asassn_data <- read.csv(asassn.csv.file,header=TRUE)
# clean obvious wild points
mag.median <- median(asassn_data$mag)
good_data <- (asassn_data$mag < (mag.median + abs(mag.cull.limit))) & (asassn_data$mag > (mag.median- abs(mag.cull.limit)))
#fit a spline
tmin <- min(asassn_data$HJD[good_data])
desmat <- asassn_data$HJD[good_data] - tmin
mars <- earth(x=desmat,y= asassn_data$mag[good_data],nk= marsOrder,pmethod= marsPMethod,penalty = marsPenalty)

# plot limits
my.y.lims <- c(max(asassn_data$mag[good_data]) + mag.margin,min(asassn_data$mag[good_data]) - mag.margin)
my.xlab <- paste("Julian Date - ",tmin)
# plot
plot(desmat,asassn_data$mag[good_data],col="darkgreen",pch=20,main="ASASSN V Data",xlab=my.xlab,ylab="V magnitude",ylim=my.y.lims)
grid(col="black")
lines(x=desmat,y=mars$fitted.values,col= "black",lwd=2)

