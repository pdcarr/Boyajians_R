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

###################################
# Confined Gaussian window
# t are the vector of sampling times (best to pass only times with valid samples)
cg.window <- function(t) {
	t.ok <- !is.na(t)
	sigma.t <- sd(t,na.rm=TRUE)
	N <- length(t[t.ok])
	w <- vector(mode="numeric")
	for (n in seq(1,N)){
#		print(n)
		w[n] <- G.x(n,sigma.t,N) - G.x(-0.5,sigma.t,N)*(G.x(n+N,sigma.t,N) + G.x(n-N,sigma.t,N))/(G.x(-0.5+N,sigma.t,N) + G.x(-0.5-N,sigma.t,N))
#		print(w[n])
	}
	return(w)
} 
##############################
# Gaussian function for Gaussian window
G.x <- function(x,sigma.t,N) {
	arg <- ((x - (N-1)/2)/(2*sigma.t))^2
#	print(arg)		# debug
	return (exp(-arg))
}