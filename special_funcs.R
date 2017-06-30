#######################################
# returns Keplerian semi-major axis in units consistent with period and mu
sma.from.period <- function(period,mu) {
	n <- 2*pi/period		# mean motion
	a <- (mu/n^2)^(1/3)
	return(a)
}