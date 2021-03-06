# R.bins
# V.bins
# B.bins
band.bias <- c(0.0,0.0,-0.380,-0.9023)
plot.title <- "R, V and B bands superimposed"
plot.pch = 20
plot.cex <- 0.5
x.pretty <- 10 # days
JD.left <- 2457900
JD.right <- 2458110
tmin = floor(JD.left/x.pretty)*x.pretty
plot.xlim <- c(0,ceil(JD.right/x.pretty)*x.pretty - tmin)
x.label <- paste("days since",tmin)
y.label <- "Magnitudes (shifted)"
# y limits
brightest.B <- min(B.bins$Magnitude) + band.bias[4]
brightest.V <- min(V.bins$Magnitude) + band.bias[3]
brightest.R <- min(R.bins$Magnitude) + band.bias[2]
dimmest.B <- max(B.bins$Magnitude) + band.bias[4]
dimmest.V <- max(V.bins$Magnitude) + band.bias[3]
dimmest.R <- max(R.bins$Magnitude) + band.bias[2]
y.limits <- c(max(dimmest.B,dimmest.V,dimmest.R),min(brightest.B,brightest.V,brightest.R))


# plot the points
quartz("superimposed")
# R band first
plot.mags = R.bins$Magnitude + band.bias[2]
plot.days <- R.bins$JD - tmin
plot(plot.days,plot.mags,
	xlim = plot.xlim,ylim=y.limits,
	main = plot.title,
	sub = "AAVSO data",
	xlab = x.label,ylab=y.label,
	type="p",col="red",pch=plot.pch,cex= plot.cex)

plot.pred <- R.bins$bin.predict + band.bias[2]
lines(plot.days, plot.pred,
	col="red",
	lwd=1)
	
# plot V
plot.mags = V.bins$Magnitude + band.bias[3]
plot.days <- V.bins$JD - tmin

points(plot.days,plot.mags,col="green",
		type="p",pch=plot.pch,cex= plot.cex)
plot.pred <- V.bins$bin.predict + band.bias[3]
lines(plot.days, plot.pred,
	col="green",
	lwd=1)


#plot B
plot.mags = B.bins$Magnitude + band.bias[4]
plot.days <- B.bins$JD - tmin

points(plot.days,plot.mags,col="blue",
		type="p",pch=plot.pch,cex= plot.cex)
		
plot.pred <- B.bins$bin.predict + band.bias[4]
lines(plot.days, plot.pred,
	col="blue",
	lwd=1)


grid(col="black")