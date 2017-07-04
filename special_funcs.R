#######################################
# returns Keplerian semi-major axis in units consistent with period and mu
sma.from.period <- function(period,mu) {
	n <- 2*pi/period		# mean motion
	a <- (mu/n^2)^(1/3)
	return(a)
}

#####################################
# returns rate of nodal precession for J2, inclination, semi-major axis, radius of primary,eccentricity, and mu
# returned units wil be consistent with mu. Inclination should be in radians
node.rate <- function(J2,sma,Rstar,incl_rad,ecc,mu) {
	semip = sma/Rstar*(1-ecc^2)	# semiparameter
	mean.motion = sqrt(mu/sma^3)
	return(-1.5*J2/semip^2*cos(incl_rad)*mean.motion)
}