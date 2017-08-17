####### load constant values for Tabby's Star in MKS units into workspace
R.Earth <- 6378e3 # radius of the Earth in meters
R.Jupiter <- 69.911e6 # radius of Jupiter
R.star <- 1.1e9	# radius if the star in meters
M.star <- 2.84e30 # mass of the star in kilograms
L.star <- 1.8e28 # luminosity of the star in Watts
G.const <- 6.674e-11 # gravitational constant Newton meters^2/kg^2
Pl.const <- 6.62607e-34 # Planck's constant in Joule seconds
mu.star <- M.star*G.const
omega.star = 2*pi/(0.88*86400) # angular rate of the star's outer atmosphere in radians per second
