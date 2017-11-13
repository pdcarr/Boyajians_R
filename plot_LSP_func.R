library(lomb)
library("MASS") # for rlm() and lqs()
#################################################
plot.LSP.period <- function(X.meas,t.times,tref=0,use.resids=TRUE,title.string=NULL,use.ofac=4) {

	desmat <- t.times - tref
	w <- cg.window(desmat)	# confined Gaussian Window
	
	if (use.resids) {
		desmat <- t.times - tref
		theFit <- rlm(X.meas ~ desmat,na.action="na.omit",psi=psi.bisquare)
		X <- as.vector(theFit$residuals*t(w))
	} else {
		X <- as.vector(X.meas*t(w))
	}
	
	mylsp <- lsp(x=theFit$residuals,times=mytimes,type="period",ofac=set.ofac,alpha=0.01)

	quartz("LSP Plot")
	plot(mylsp$scanned,mylsp$power,xlab="Period (days)",ylab="power",col="darkgreen",type="l",main=title.string)

}