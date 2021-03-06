dips <- NULL
############# Inputs #############
klic.file <- "data/kplr008462852-2009350155506_llc_t1.txt"
klic.file <- "data/kplr008462852-2013131215648_llc_t1.txt"
#klic.file <- "data/kplr008462852-2012179063303_llc_t1.txt"
#klic.file <- "data/kplr008462852-2009131105131_llc_t1.txt" # QUARTER 0
#klic.file <- "data/kplr008462852-2013098041711_llc_t1.txt"

plot.title <- "Kepler PDCSAP data"
earliest.day <- 0
#earliest.day <- 1470
#earliest.day <- 1200
last.day <- 1600
#last.day <- 1560
#last.day <- 1220
NA.string <- "NA"
NA.string <- "NULL"
plot_both = FALSE
pad = 100

#dips <- data.frame(dstart=1204,dstop =1209,stringsAsFactors=FALSE)
dips <- data.frame(dstart=1505,dstop=1525,stringsAsFactors=FALSE)
dips <-  rbind(dips,c(1530,1550))

###################
klightcurve <- read.csv(file=klic.file,header=FALSE,check.names=TRUE,na.strings=NA.string,stringsAsFactors=FALSE)

colnames(klightcurve) <-  c("Time","TIMECORR","CADENCENO","SAP_flux","SAP_flux_err","SPA_BKG","SAP_BKG_ERR","PDCSAP_flux","PDCSAP_flux_err","SAP_quality","PSF_centr1","PSF_centr1_err",
	"PSF_CENTR2","PSF_CENTR2_ERR","MOM_CENTR1","MOM_CENTR1_ERR","MOM_CENTR2","MOM_CENTR2_ERR","POS_CORR1","POS_CORR2")
	
# apply time bounds

in.bounds <- klightcurve$Time <= last.day & klightcurve$Time >= earliest.day

in.bounds <- in.bounds & klightcurve$SAP_quality != 2176	 & klightcurve$SAP_quality != 8192	# cosmic ray or other fault
	
#
myuppery <- max(klightcurve$PDCSAP_flux[in.bounds],na.rm=TRUE)
mylowery <- min(klightcurve$PDCSAP_flux[in.bounds],na.rm=TRUE)
if(plot_both) {
	myuppery <- max(myuppery,max(klightcurve$SAP_flux[in.bounds],na.rm=TRUE),na.rm=TRUE)
	mylowery <- min(mylowery,min(klightcurve$SAP_flux[in.bounds],na.rm=TRUE),na.rm=TRUE)
}

myylims <- c(mylowery - pad,myuppery + pad)

plot(klightcurve$Time[in.bounds],klightcurve$PDCSAP_flux[in.bounds],type="l",col="darkblue",xlab="JD - 2454833",ylab="flux",ylim=myylims,main=plot.title)
if (plot_both) {
	lines(klightcurve$Time[in.bounds], klightcurve$SAP_flux[in.bounds],col="red")
}
grid(col="black")