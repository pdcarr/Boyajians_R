rm("sgb")
if(exists('binCurve') & exists('tmin')) {
	sgb <- binCurve$Band == "SG"
	if(sum(sgb) <= 0) { 
		cat("\nno SG data in binCurve\n")
	} else {		
		
		SO.time <- binCurve$JD[sgb] - 2400000.5 - tmin
		points(x=SO.time,
			y=binCurve$Magnitude[sgb],
			col="aquamarine3",
			pch=1,cex=2.0)
	}
}