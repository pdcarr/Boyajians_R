blue <- lightcurve$Band == "B"
x <- unique(lightcurve$Observer_Code[blue])
print("Code,N_obs,airmass_cnt,comp_1,comp_2")
for (observer in x) {
	mine <- lightcurve$Observer_Code == observer & blue
	nobs <- length(lightcurve$JD[mine])
	good_air_cnt <- length(lightcurve$JD[!is.na(lightcurve$Airmass) & mine])
	cstars.1 <- unique(lightcurve$Comp_Star_1[mine])
	cstars.2 <- unique(lightcurve$Comp_Star_1[mine])
	obs.string <- paste(observer,nobs,good_air_cnt,cstars.1,cstars.2)
	print(obs.string)
}