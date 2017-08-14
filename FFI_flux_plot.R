
library("smooth") # for smoothing
library("earth") # for MARS
library("crayon") # to add color to text
# read in data file from Montet and Simon paper (cleaned up a bit and converted to .csv)
MandSFFI <- read.csv(file="data/Montet_SIMON_FFI_table.csv",header=TRUE)
# perform a linear spline fit increase nk and decrease penalty to see more features, or set set pmethod = "none"
thisFit <- earth(x= MandSFFI$Time,y= MandSFFI$Normalized.Flux,nk= 7,pmethod= "backward",penalty = 3)
# plot this stuff
plotTitle <- paste(c("Normalized Flux Data from Kepler FFI Images","http://iopscience.iop.org/article/10.3847/2041-8205/830/2/L39"),collapse="\n")

plot(x=MandSFFI$Time,y=MandSFFI$Normalized.Flux,col="red",xlab="JD - 2454833",ylab="Normalized Flux",main=plotTitle,cex.main=0.7,pch=3)
lines(x= MandSFFI$Time,y=thisFit$fitted.values,col="black",lwd=2)
grid(col="black")