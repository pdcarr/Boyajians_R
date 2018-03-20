dip.mask <- data.frame(dipname="Elsie",JD.begins = 2457890,JD.ends=2457896,stringsAsFactors=FALSE)
dip.mask <- rbind(dip.mask,c("Celeste",2457915,2457927),
						   c("Skara Brae",2457968,2457984),
						   c("Angkor",2457994,2458011),
						   c("Ozymandias",2458194,2458197))

# adding mask of FWAIN, as observed by Bruce Gary.
#dip.mask <-rbind(dip.mask,c("FWAIN",2458051,2458057))

# primarily observed by Gary Johnson, the "December Surprise"
#dip.mask <- rbind(dip.mask,c("December_Surprise",2458096,2458097))
