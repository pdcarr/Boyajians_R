###########################################################
AirMass <- function(JD,locObs,starLoc) {
# JD is the vector of Julian dates
# locObs is the decimal location = c(lat,long) of the observatory
# starLoc is a vector of declination degrees, minutes, seconds and right ascension in h,m,s of the star
# calculates Airmass from the time, the observer's location, and the declination of the star
# uses astroFns package
	library(astroFns)
	
	# modified Julian Date
	mJD <- JD - 2400000.5
	
	# calculate UT times
	utStrings <- dmjd2ut(mJD,tz="UTC")
	
	#break out the elements of the time as a POSIXlt class object
	ltTimes <- as.POSIXlt(utStrings)
	# calculate hour angles at the observatory
	myHAs <- ut2ha(yr=ltTimes$year,mo=ltTimes$mon + 1,dy=ltTimes$mday,hr=ltTimes$hour,mi=ltTimes$min,se=ltTimes$sec,ra.sou =starLoc[2],lon.obs=rad2hms(locObs[2]*pi/180))
#	print(myHAs)
	#calculate elevation angles from hour angles and observatory latitude
	myEls = elev(dec.sou=starLoc[1],ha=myHAs,lat.obs=rad2dms(locObs[1]*pi/180))
#	print(myEls)
	return(1/abs(sin(myEls*pi/180)))
}

########################################################
ReverseMagnitude <- function(relMag) {
	# relMag is a difference in astronomical magnitudes
	return(10^(-relMag/2.5))
}
