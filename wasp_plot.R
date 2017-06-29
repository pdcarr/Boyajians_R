# read in, fit and plot the wasp data
library(MASS)

wasp_data <- read.csv("wasp.csv",header=FALSE)

plot(wasp_data$V2,wasp_data$V3,col="red",pch=20,cex=0.7,ylim=c(max(wasp_data$V3,na.rm=TRUE),min(wasp_data$V3,na.rm=TRUE)),xlab="Julian Date",ylab="Magnitude",main="WASP Data")

myfit = rlm(wasp_data$V3 ~ wasp_data$V2)
lines(wasp_data$V2,myfit$fitted.values,col="black")

grid(col="black")
