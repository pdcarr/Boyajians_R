# fill in some missing airmass values
lasCruces <- c(32.31994,-106.763654) # center of Las Cruces, NM in decimal degrees latitude, longitude.
Leominster <- c(52.226529,-2.741) # approximate location of PXR
Warrington <- c(53.39,-2.59695) # approx. location of JSJA
Erd <- c(47.39,18.90454) # approx. location of MATA
Haan <- c(51.192,7.007) # approx. location of WROC
LosAngeles <- c(34.05,-118.25) # approx. location of FJQ
Potomac <- c(38.02,-77.21) # approx location of PALE
Coventry <- c(52.4,-1.52) # approximate location of KTHC
Ponferrada <- c(42.55,-6.598) # approximate location of PJVA
Cotopaxi <- c(38.38,-105.69) # approximate location of HJW
Minneapolis <- c(44.9778,-93.265)
Auburn <- c(41.367,-85.06) # approximate location for SDB
Sorrento <- c(40.63,14.38) # approximate location ofr RNL
Budapest <- c(47.47,19.04) # approximate location for TIA
ValRacine <- c(45.48,-71.08) # approximate location of DJED
France <- c(46.3,2.2)  # very rough location of SFLB
tabbysLoc <- c("+44d 27m 24.61s","20h 06m 15.457s") # declination and right ascension of the star.
missingAirmass <- c("JM","PXR","JSJA","MATA","WROC","FJQ","PALE","KTHC","PJVA","HJW","DRZ","SDB","RNL","TIA","DJED","SFLB")	# observer code
missingAMLocs = rbind(lasCruces,Leominster,Warrington,Erd,Haan,LosAngeles,Potomac,
						Coventry,Ponferrada,Cotopaxi,Minneapolis,Auburn,Sorrento,Budapest,
						ValRacine,France)
missingAM.data <- data.frame(lax.observers = missingAirmass,missing.lats = missingAMLocs[,1],missing.lons=missingAMLocs[,2],stringsAsFactors = FALSE)