library("MASS") # for rlm() and lqs()
library("Hmisc")
source("data_funcs.R")
rm("lco.data")
########## setup paramters
epsilon.t = 1/86400 # small time interval to start binning
##########
trial.bins <- 200 # of bins to try in kmeans
#t.margin <- 1 # days
############## exclude any observatories or cameras (case sensitive)?
exclude.codes <- c("TFN")
########### set up the file(s)
#LCO_file <- "data/Measurements_subset_I_TFN_KB25.csv"

LCO.bands <- c("I","B")
data.directory <- "../LCO_GS/"
LCO.suffix <-"txt$"
LCO.prefix <- "^Measurements_subset_"

# figure out what files are available and load them in
all.the.files <- dir(data.directory)
these.files <- grepl(pattern=LCO.suffix,x=all.the.files,perl=TRUE)
lco.data <- data.frame(Num = numeric(),label=character(),MJD=numeric(),rel_flux_T1_n=numeric(),band=character()) # create data fram for all the data
file.data <- data.frame(Num = numeric(),label=character(),MJD=numeric(),rel_flux_T1_n=numeric()) # create scratch data fram for the data read in

for (this.band in LCO.bands) {
  band.pattern <- paste(LCO.prefix,this.band,".*",LCO.suffix,sep = "")
  band.files <- grepl(pattern=band.pattern,x=all.the.files[these.files],perl=TRUE)
  for (this.file in all.the.files[these.files][band.files]) {
    this.file <- paste(data.directory,this.file,sep="")
    file.data <- read.delim(file=this.file,header=TRUE,col.names=names(file.data),stringsAsFactors=FALSE)
    nrec = length(file.data[,1])
    band.column <- data.frame(band=rep(this.band,nrec))
    lco.data <- rbind(lco.data,cbind(file.data,band.column))
  }
}

#LCO.files <- rbind(LCO.files,))

#########
plot.col = "darkviolet"
plot.sym = 20


for(myband in unique(lco.data$band)) {
  quartz("LCO data")
	these.obs <- lco.data$band == myband
	my.title <- paste("LCO ",myband,"Data - ",sum(these.obs),"Measurements")
	plot(x=lco.data$MJD[these.obs],
	     y=lco.data$rel_flux_T1_n[these.obs],
	     type="p",
	    col=plot.col,pch=plot.sym,
		  xlab="JD - 2400000",ylab="normalized flux",
		  main=my.title)
	grid(col="black")
# pick bins using K Means (dirt simple algortithm)
	bin.data <- kmeans.time.series(times=lco.data$MJD[these.obs],initial.clusters.num=trial.bins,min.population=2,delta.mean=0.01,max.iterations=12)
	# our.bins <- bin.data[1]
	# mem <- bin.data[2]
	bin.numbers <- seq(1,length(bin.data$bins),1)
	bin.flux.means <- sapply(bin.numbers,bin.fluxes,bin.data$membership,lco.data$rel_flux_T1_n[these.obs])
	quartz("binned plot of LCO data")
	plot(x=lco.data$MJD[these.obs],
	     y=lco.data$rel_flux_T1_n[these.obs],
	     type="p",
	     col=plot.col,pch=plot.sym,
	     xlab="JD - 2400000",ylab="normalized flux",
	     main=my.title)
	grid(col="black")
	
	
#	hold on
}
