#### Airmass() ####################################################### 
AirMass <- function(JD,locObs,starLoc) {
# JD is the vector of Julian dates
# locObs is the decimal location = c(lat,long) of the observatory
# starLoc is a vector of declination degrees, minutes, seconds and right ascension in h,m,s of the star
# calculates Airmass from the time, the observer's location, and the declination of the star
# uses astroFns package
	library(astroFns)
#	print(locObs)
	# modified Julian Date
	mJD <- JD - 2400000.5
	
	# calculate UT times
	utStrings <- dmjd2ut(mJD,tz="UTC")
#	print(tail(utStrings))
	
	#break out the elements of the time as a POSIXlt class object
	ltTimes <- as.POSIXlt(utStrings)
	# calculate hour angles at the observatory
	myHAs <- ut2ha(yr=ltTimes$year +1900,mo=ltTimes$mon + 1,dy=ltTimes$mday,hr=ltTimes$hour,mi=ltTimes$min,se=ltTimes$sec,
					ra.sou =starLoc[2],lon.obs=rad2dms(locObs[2]*pi/180))
# 	print(tail(ltTimes))
#	print(tail(myHAs))
	#calculate elevation angles from hour angles and observatory latitude
	myEls = elev(dec.sou=starLoc[1],ha=myHAs,lat.obs=rad2dms(locObs[1]*pi/180))
#	print(myEls)
	return(1/abs(sin(myEls*pi/180)))
}

#### ReverseMagnitude() ####################################################
ReverseMagnitude <- function(relMag) {
	# relMag is a difference in astronomical magnitudes
	return(10^(-relMag/2.5))
}

#### HillRadius() ######################################################
# returns approximate Hill radius in units consistent with a for mass of secondary (m) and Primary (M)
HillRadius <- function(m=1,M=1000,a,e=0) {
	return(a*(1-e)*(m/(3*M))^(1/3))
}