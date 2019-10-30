	cat("\n Welcome to Boyajian's R!\n")
	cat("\nsetting 	options(digits=12)\n")
	options(digits=12)
	if(exists("lightcurve") & is.data.frame(lightcurve)) {
	  notas <- lightcurve$Observer_Code != "ASASSN"
	  print(aavsoJDbounds(lightcurve[notas,]))
	} else {
	  cat("\nNo light curve data frame found\n")
	}
