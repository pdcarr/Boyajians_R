options(digits=20)
library(lattice)

source.id <- 2081900940499099136


#data.file.name <- "data/GaiaSource_000-000-000.csv"
#data.file.name <- "data/GaiaSource_000-004-066.csv"
#data.file.name <- "data/GaiaSource_000-003-200.csv"
#data.file.name <- "data/GaiaSource_000-003-184.csv"
data.file.name <- "data/GaiaSource_000-003-178.csv"


#data.file.name <- "data/GaiaSource_000-020-001.csv"
#data.file.name <- "data/GaiaSource_000-020-081.csv"
#data.file.name <- "data/GaiaSource_000-020-110.csv"
cat("reading file\n")
gaia <- read.csv(file = data.file.name, header = TRUE, stringsAsFactors = FALSE)
cat("plotting\n")
plot(gaia$ra, gaia$dec, pch = 20, type = "p", cex = 0.1, xlab = "RA (deg)", ylab = "dec(deg)", col = "red")
grid(col="black")

sid.min <- min(gaia$source_id)

sid.max <- max(gaia$source_id)

if(source.id >= sid.min & source.id <= sid.max) {
	cat("\n target source id is in range")
	my.records <- gaia$source_id == source.id
	hits <- length(my.records)
	print(paste(hits,"records match source id ",source.id))
	} else {
	cat("\n source id not in range")
}

