###################
weightedBins <- TRUE # set to TRUE to weight lower uncertainty bins more and exclude bins in dips

####### MARS
plotMARS <-  FALSE # set to TRUE to try a MARS fit instead of lm() or rlm()
# fit parameters
marsOrder <- 21 # maximum number of knots
marsPenalty <- 5 # set to 0 to avoid penalizing knots in pruning pass
#marsPMethod <- "none" # set to "none" to avoid pruning
marsPMethod <- "exhaustive" # set to "none" to avoid pruning
mars.thresh <- 0.00002 # threshold parameter for earth().
mars.minspan <- 5 # minimum number of observations between knots
splineRaw <-  FALSE # do the spline on the raw lightcurve, not binned.
mask.Dips <- TRUE

############# smooth.spline
perform.smooth <- TRUE
smooth.n.knots <- 9
df.penalty <- 1
do.CV <- TRUE
smoove.color <- "slategray"
