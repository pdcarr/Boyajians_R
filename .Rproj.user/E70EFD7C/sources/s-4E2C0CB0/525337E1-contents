###################
weightedBins <- TRUE # set to TRUE to weight lower uncertainty bins more and exclude bins in dips

####### MARS
plotMARS <-  FALSE # set to TRUE to try a MARS fit instead of lm() or rlm()
# fit parameters
marsOrder <- 7 # maximum number of knots
marsPenalty <- 6 # set to 0 to avoid penalizing knots in pruning pass
#marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "exhaustive" # set to "none" to avoid pruning
mars.thresh <- 0.00002 # threshold parameter for earth().
mars.minspan <- 2 # minimum number of observations between knots
splineRaw <-  FALSE # do the spline on the raw lightcurve, not binned.
mask.Dips <- TRUE
############# smooth.spline
perform.smooth <- TRUE # TRUE if you want to plot a smooth cubic spline
smooth.deriv <- TRUE	 # TRUE if you want to plot the first derivative of the smooth spline as well.
deriv.margin = 1.02
smooth.n.knots <- 6 # number of knots, including the two endpoints
df.penalty <- 1
do.CV <- TRUE	# do cross validation
smoove.color <- "slategrey" # color of the line drawn on the magnitude plot
