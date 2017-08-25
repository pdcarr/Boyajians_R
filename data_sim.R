library(lomb)
### simulate a lot of noise and a little sinusoid
N.pts <- 1000
sin.amp = 0.001
sin.period = 10 # keep well under Nyquist
sin.phase <- 0
ofac.par <- 8

################
the.noise <- arima.sim(n=N.pts,model=list(c(0,0,0)))
t <- seq(1,N.pts)
X <- vector()
for (index in t) {
 X[index] <- the.noise[index] + sin((2*pi/sin.period)*t[index] + sin.phase)
}

quartz("simulated data")
plot(t,X,main="simulated data",xlab="index",ylab="signal",type="l")
grid(col="black")

#


quartz()
sim.lsp <- lsp(x=X,times=t,type="period",ofac=ofac.par,alpha=0.01)
quartz("Simulated LSP")
title <- paste("simulated LSP with sin magnitude =",sin.amp)
plot(sim.lsp$scanned,sim.lsp$power,xlab="Period",ylab="normalized power",col="red",type="l",main=title,xlim=c(0,2*sin.period),log="y")
#par(xlog=TRUE)
grid(col="black")
#
window.string <- "Confined Gaussian Window"
quartz("simulated with window")
w <- cg.window(t)
sim.w.lsp <- lsp(x=as.vector(X*t(w)),times=t,type="period",ofac=ofac.par,alpha=0.01)
title <- paste("simulated LSP with sin magnitude =",sin.amp,"with",window.string)
plot(sim.w.lsp$scanned, sim.w.lsp$power,xlab="Period",ylab="normalized power",col="red",type="l",main=title,xlim=c(0,2*sin.period),log="y")
grid(col="black")



