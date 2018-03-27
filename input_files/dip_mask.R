dip.mask <- data.frame(dipname="Elsie",JD.begins = 2457890,JD.ends=2457896,stringsAsFactors=FALSE)
dip.mask <- rbind(dip.mask,c("Celeste",2457915,2457927),
						   c("Skara Brae",2457968,2457984),
						   c("Angkor",2457994,2458011),
						   c("Ozymandias",2458191.9,2458206))

# adding mask of FWAIN, as observed by Bruce Gary.
#dip.mask <-rbind(dip.mask,c("FWAIN",2458051,2458057))

# primarily observed by Bruce Gary, the "December Surprise"
dip.mask <- rbind(dip.mask,c("December_Surprise",2458094,2458097))
