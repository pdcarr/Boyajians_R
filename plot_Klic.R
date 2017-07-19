
klic.file <- "data/kplr008462852-2009350155506_llc_t1.txt"
plot_both = FALSE
pad = 100
###################
klightcurve <- read.csv(file=klic.file,header=FALSE,check.names=TRUE,na.strings="NA")

colnames(klightcurve) <-  c("Time","TIMECORR","CADENCENO","SAP_flux","SAP_flux_err","SPA_BKG","SAP_BKG_ERR","PDCSAP_flux","PDCSAP_flux_err","SAP_quality","PSF_centr1","PSF_centr1_err",
	"PSF_CENTR2","PSF_CENTR2_ERR","MOM_CENTR1","MOM_CENTR1_ERR","MOM_CENTR2","MOM_CENTR2_ERR","POS_CORR1","POS_CORR2")
	
#
myuppery <- max(klightcurve$PDCSAP_flux,na.rm=TRUE)
mylowery <- min(klightcurve$PDCSAP_flux,na.rm=TRUE)
if(plot_both) {
	myuppery <- max(myuppery,max(klightcurve$SAP_flux,na.rm=TRUE),na.rm=TRUE)
	mylowery <- min(mylowery,min(klightcurve$SAP_flux,na.rm=TRUE),na.rm=TRUE)
}

myylims <- c(mylowery - pad,myuppery + pad)

plot(klightcurve$Time,klightcurve$PDCSAP_flux,type="l",col="darkblue",xlab="JD - 2454833",ylab="flux",ylim=myylims)
if (plot_both) {
	lines(klightcurve$Time, klightcurve$SAP_flux,col="red")
}
grid(col="black")