library("MASS") # for rlm() and lqs()
library("Hmisc")
library("pracma")
### this script uses the data_funcs library
source("data_funcs.R")
######## remove old stuff to avoid confusion
rm("lco.data")
rm("bin.data")
rm("resistFit")
options(digits=12)
########## setup paramters
regression.methods <- c("None",
                        "resistant",
                        "robust")
method.2.use <- 3
##########
trial.bins <- 200 # of bins to try in kmeans
#t.margin <- 1 # days
############## exclude any observatories or cameras (case sensitive)?
exclude.codes <- c("TFN") # not yet implemented

############## filter band and plot setup
plot.raw <- FALSE # TRUE if you want to plot data for each band right frm the files with no binning
# use this data frame to set up colors and symbols for all the band you expect to be in the files.
LCO.bands <- data.frame(band.codes=c("R","B","I"),band.colors=c("red","blue","darkviolet"),band.symbol=c(20,20,20),stringsAsFactors=FALSE)
# LCO.bands <- data.frame(band.codes=c("I"),band.colors=c("darkviolet"),band.symbol=c(20),stringsAsFactors=FALSE) # I band only
# LCO.bands <- data.frame(band.codes=c("B"),band.colors=c("blue"),band.symbol=c(20),stringsAsFactors=FALSE) # B band only
############### set up plot time bounds
earliest_MJD <- NA
# earliest_MJD <- 58770
pretty.days <- 5

########### set up the search for file(s)

data.directory <- "../LCO_GS/"
LCO.suffix <- "(\\h(TFN|ELP|OGG)\\h|_).*txt$"
LCO.prefix <- "(^Measurements_subset_|^)"


# figure out what files are available and load them in
all.the.files <- dir(data.directory)
these.files <- grepl(pattern=LCO.suffix,x=all.the.files,perl=TRUE) # logical vector of matching data files
lco.data <- data.frame(Num = numeric(),label=character(),MJD=numeric(),rel_flux_T1_n=numeric(),band=character()) # create data fram for all the data
file.data <- data.frame(Num = numeric(),label=character(),MJD=numeric(),rel_flux_T1_n=numeric()) # create scratch data frame for the data read in
# loop over bands, then the files matching each band
for (this.band in LCO.bands$band.codes) {
  band.pattern <- paste(LCO.prefix,this.band,LCO.suffix,sep = "") # put the search pattern together for the band
  band.files <- grepl(pattern=band.pattern,x=all.the.files[these.files],perl=TRUE) # all the matching files for the band in question
  for (this.file in all.the.files[these.files][band.files]) {
    this.file <- paste(data.directory,this.file,sep="")
    file.data <- read.delim(file=this.file,header=TRUE,col.names=names(file.data),stringsAsFactors=FALSE)
    nrec = length(file.data[,1]) # number of records in the file
    band.column <- data.frame(band=rep(this.band,nrec)) # construct a column to be bound to the data for that band
    lco.data <- rbind(lco.data,cbind(file.data,band.column)) # glue the band column unto the data from the file.
  }
}

t.min.LCO <- min(lco.data$MJD,na.rm=TRUE) # the earliest time in the data frame

#########

# loop over the bands
for(myband in unique(lco.data$band)) {
  band.index <- match(myband,LCO.bands$band.codes)
# browser()
    if(!is.na(band.index)) {
      plot.col <- LCO.bands$band.colors[band.index] # the plot color for this band
      plot.sym <- LCO.bands$band.symbol[band.index] # the plot symbol for this band
  } else {
      print(paste("Band symbol",myband,"not recognized"))
      next()
  }

  these.obs <- lco.data$band == myband
  if(plot.raw) {
    quartz("LCO data")
  	my.title <- paste("LCO ",myband,"Data - ",sum(these.obs),"Measurements")
# x axis time limits
  	latest.MJD <- pretty.days*ceil(max(lco.data$MJD[these.obs])/pretty.days)
  	if(is.na(earliest_MJD)){
  	 plot.start.MJD <- pretty.days*(floor(min(lco.data$rel_flux_T1_n[these.obs])/pretty.days))
  	} else {
  	  plot.start.MJD <- earliest_MJD
  	}
  	plot.x.lims <- c(plot.start.MJD,latest.MJD)
# browser() 
  	plot(x=lco.data$MJD[these.obs],
  	     y=lco.data$rel_flux_T1_n[these.obs],
  	     type="p",
  	    col=plot.col,pch=plot.sym,
  	    xlim=plot.x.lims,
  		  xlab="JD - 2400000",ylab="normalized flux",
  		  main=my.title)
  	grid(col="black")
  }
########## bin and plot binned data for the color in question (myband)  
# pick bins using K Means with pruning (dirt simple algorithm)
  # browser()
	bin.data <- kmeans.time.series(times=lco.data$MJD[these.obs],initial.clusters.num=trial.bins,min.population=2,delta.mean=0.01,max.iterations=12)
	# our.bins <- bin.data[1]
	# mem <- bin.data[2]
	bin.numbers <- seq(1,length(bin.data$bins),1) # a sequence from 1 to the number of bins
	bin.flux.means <- sapply(bin.numbers,bin.fluxes,bin.data$membership,lco.data$rel_flux_T1_n[these.obs]) # mean flux for each bin
	tmin <- min(bin.data$bins,na.rm=TRUE) # earliest bin time
	
	# calculate a regression
	desmat <- bin.data$bins - tmin  # subtract off the earliest time to make the coeeficients easier to interpret
	# use resistant linera regression 
#	resistFit <- lqs(formula = bin.flux.means ~ desmat)
	quartz("binned plot of LCO data")
	bin.title <- paste(c(paste("binned measurements",myband,"Band"),paste("regression=",regression.methods[method.2.use])),collapse="\n")
# set X axis time limits
	latest.MJD <- pretty.days*ceil(max(bin.data$bins)/pretty.days)
	if(is.na(earliest_MJD)){
	  plot.start.MJD <- pretty.days*(floor(min(bin.data$bins)/pretty.days))
	} else {
	  plot.start.MJD <- earliest_MJD
	}
	plot.x.lims <- c(plot.start.MJD,latest.MJD)
# plot!	
	plot(x=bin.data$bins,
	     y=bin.flux.means,
	     type="p",
	     col=plot.col,pch=plot.sym,
	     xlim=plot.x.lims,
	     xlab="Bin center MJD",ylab="mean normalized flux",
	     main=bin.title)
	grid(col="black")
	
###### pick and execute the regression method
	
	if (regression.methods[method.2.use] == "None") {
	  cat("\nNo regression method selected\n")
	} else if(regression.methods[method.2.use] =="resistant") {

    	resistFit <- lqs(formula = bin.flux.means ~ desmat);
    	lines(x=resistFit$model$desmat + tmin,
    	      y=resistFit$fitted.values,
    	      col="darkgrey",lwd=2);
    	cat("color: ",myband,"\ncoefficients: \n");
    	print(coefficients(resistFit));
    	cat("\n") 
    	
	} else if(regression.methods[method.2.use] == "robust") {
	    robustFit <- rlm(bin.flux.means ~ desmat,na.action="na.omit",psi=psi.bisquare)
	    lines(x=robustFit$model$desmat + tmin,
	          y=robustFit$fitted.values,
	          col="darkgrey",lwd=2);
	    cat("color: ",myband,"\ncoefficients: \n");
	    print(coefficients(robustFit));
	    cat("\n") 
	    SE.slope <- slope.SE(robustFit$model$desmat,robustFit$residuals)
	    print(paste("Slope standard error: ",SE.slope))
	    
	}else {
    cat("\n no regression method specified\n")
  }

#	hold on
}
