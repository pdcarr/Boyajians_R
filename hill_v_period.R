
source("input_files/facts_MKS.R")
source("astro_funcs.R")
mass.planet <- 10 # Jupiter masses
m.star.jup <- M.star/M.Jupiter
period <- seq(1,50,0.5) # period in days
ecc <- 0
hr <- vector()
for (P in period) {
	mm <- 2*pi/(P* seconds.ms.day)
	sma <- (mu.star/mm^2)^(1/3)
	hr <- rbind(hr,HillRadius(10, m.star.jup,sma,ecc))
}
title <- paste("Mass of planet =",mass.planet,"MJ")
plot(period,hr/(R.star),xlab="Period (days)",ylab="Hill Radius (stellar radii)",type="l",main=title)
grid(col="black")