library("MASS") # for rlm() and lqs()
library("Hmisc")
rm("lco.data")
rm("band.vec")
########### set up the file(s)
#LCO_file <- "data/Measurements_subset_I_TFN_KB25.csv"
rm("LCO.files")
LCO.files <- data.frame(lco.file="data/Measurements_subset_I_TFN_KB25.csv",band="I",stringsAsFactors=FALSE)
LCO.files <- rbind(LCO.files, c("data/Measurements_subset_I_TFN_KB23.csv","I"),c("data/Measurements_subset_I_OGG_KB82.csv","I"))
##########
bin.width <- 1 # days
t.margin <- 1 # days
#########
plot.col = "darkviolet"
plot.sym = 20

# open the file(s)
got.good.file <- FALSE
f.index <- 1:nrow(LCO.files)
lco.data <- data.frame(Num = numeric(),label=character(),JD.2400000=numeric(),rel_flux_T1_n=numeric(),band=character())
for(this.file in f.index) {
	lco.scratch <- NA
	lco.scratch <- read.csv(file= as.character(LCO.files$lco.file[this.file]),stringsAsFactors=FALSE,header = TRUE)
	if(!is.na(lco.scratch)) {
		if(!got.good.file) {
			got.good.file <- TRUE
			band.vec <- vector(mode="character")
			band.vec[1:nrow(lco.scratch)] <- LCO.files[this.file,"band"]
			lco.scratch <- cbind(lco.scratch,band.vec)
			lco.data <- lco.scratch
		} else {
			band.vec <- vector(mode="character")
			band.vec[1:nrow(lco.scratch)] <- LCO.files[this.file,"band"]
			lco.scratch <- cbind(lco.scratch,band.vec)
			lco.data <- rbind(lco.data,lco.scratch)
		}
	} else {print("problem reading file")}
}
quartz("LCO data")

for(myband in unique(lco.data$band)) {
	these.obs <- lco.data$band == myband
	plot(x=lco.data$JD.2400000[these.obs],y=lco.data$rel_flux_T1_n[these.obs],
		col=plot.col,pch=plot.sym,
		xlab="JD - 2400000",ylab="normalized flux")
	grid(col="black")
#	hold on
}
t.min <- min(lco.data$JD.2400000,na.rm=TRUE)