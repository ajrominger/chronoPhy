## script to run exp la and mu rate model over different par for 400 gen
setwd('~/Dropbox/Research/chronoPhy')
source('chronoPhy_nonLinear3.R')

library(parallel)


# curve(la0*exp(-bla*x), from=0, to=5, col='blue', lwd=3,
      # xlab='Time', ylab='Per species rates')
# curve(mu0*exp(-bmu*x), add=TRUE, col='red', lwd=3)
# abline(h=nu0, col='gray35', lwd=3)
# abline(h=ga0, col='blue', lwd=3, lty=2)

# x <- evolve(LaFun=RhoExpR, mu0=mu0, bmu=bmu, nu0=1, ga0=1, nGen=400, rho0=la0, b=bla)
# plot(evol2phylo(x))
# sum(rowSums(x[[1]], na.rm=TRUE) == 0)
# colSums(x[[1]], na.rm=TRUE)


## set rate parameters
la0 <- 1.4
bla <- 0.75
mu0 <- la0*exp(-bla*5)
bmu <- -1.1*bla
nu0 <- 1 #la0*0.25
ga0 <- 0.8 #nu0*0.9

la.col <- hsv(0.7, 0.8, 0.7)
mu.col <- hsv(0.1, 0.4)
nu.col <- hsv(0.5, 0.5, 0.7)

## see how they look
pdf(file='poster/fig_expLaMuRates.pdf', width=4, height=4)
par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
curve(la0*exp(-bla*x), from=0, to=5, col=la.col, lwd=3,
      xlab='', ylab='',
      axes=FALSE, frame.plot=TRUE)
curve(mu0*exp(-bmu*x), add=TRUE, col=mu.col, lwd=3)
abline(h=nu0, col=nu.col, lwd=3)
abline(h=ga0, col=la.col, lwd=3, lty=2)
dev.off()

pdf(file='poster/fig_expLaMuRates_skinny.pdf', width=6, height=5)
par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
curve(la0*exp(-bla*x), from=0, to=5, col=la.col, lwd=3,
      xlab='', ylab='',
      axes=FALSE, frame.plot=TRUE)
curve(mu0*exp(-bmu*x), add=TRUE, col=mu.col, lwd=3)
abline(h=nu0, col=nu.col, lwd=3)
abline(h=ga0, col=la.col, lwd=3, lty=2)
dev.off()


## replicate trees
nrep <- 500

noReturn <- mclapply(1:nrep, mc.cores=20, FUN=function(i) {
	x <- evolve(LaFun=RhoExpR, mu0=mu0, bmu=bmu, nu0=1, ga0=1, nGen=400, rho0=la0, b=bla)
	
	if(!is.null(nrow(x[[1]]))) {
		if(any(colSums(x[[1]], na.rm=TRUE) != 0)) {
			save(x, file=sprintf('expLaMuOut/simOut_constantRate_%s.RData', i))
		}
	}
})