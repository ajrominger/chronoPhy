## script to run exp la and mu rate model over different par for 400 gen
setwd('~/Dropbox/Research/chronoPhy')
source('chronoPhy_nonLinear3.R')

library(parallel)


## set rate parameters
la0 <- 1.4
bla <- 0.3
mu0 <- 0.1*la0
bmu <- -0.45
nu0 <- la0*0.4
ga0 <- nu0*0.8

## see how they look
la.col <- hsv(0.7, 0.8, 0.7)
mu.col <- hsv(0.1, 0.4)
nu.col <- hsv(0.5, 0.5, 0.7)
pdf(file='poster/fig_quadLaRates.pdf', width=4, height=4)
par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
curve(la0*exp(-bla*(x-1.5)^2), from=0, to=5, col=la.col, lwd=3,
      xlab='', ylab='',
      axes=FALSE, frame.plot=TRUE)
curve(mu0*exp(-bmu*x), add=TRUE, col=mu.col, lwd=3)
abline(h=nu0, col=nu.col, lwd=3)
abline(h=ga0, col=la.col, lwd=3, lty=2)
dev.off()

## replicate trees
nrep <- 500

noReturn <- mclapply(1:nrep, mc.cores=20, FUN=function(i) {
	x <- evolve(LaFun=RhoQuadR, mu0=mu0, bmu=bmu, nu0=nu0, ga0=ga0, nGen=400, rho0=la0, b=bla)
	
	if(!is.null(nrow(x[[1]]))) {
		if(any(colSums(x[[1]], na.rm=TRUE) != 0)) {
			save(x, file=sprintf('quadLaOut/simOut_constantRate_%s.RData', i))
		}
	}
})