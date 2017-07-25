######################## draw vertical lines at various dates
# Arguments: 
#	jdLine: a list of x-axis values on which to draw these dashed lines
# 	jdLineText: a list of text labels of same length as jdLine
#	plot.y.Lims:  the lower and upper y Limits for the plot.
#	jdLineColor: the color you want for the line(s)
#	jdcex: size reduction factor for text
#	jdlty: line type
verticalDateLines <- function(jdLine,jdLineText,plot.y.Lims,jdLineColor="red",jdcex=0.5,jdlty = "dashed") {

	# draw a vertical line for a date of interest
	for(jd.index in 1:length(jdLine)) {
		lines(x=c(jdLine[jd.index],jdLine[jd.index]),y=plot.y.Lims,col=jdLineColor,lwd=1,lty="dashed")
		text(x= jdLine[jd.index],y=plot.y.Lims[2],labels=jdLineText[jd.index],pos=3,cex=jdcex)
	}
}

#########################
