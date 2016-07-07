## script to run constant rate model over different par for 1000 gen

setwd('~/Dropbox/Research/chronoPhy')
source('chronoPhy.R')

library(parallel)

## set up simulation par
nGen <- 1000
nRep <- 100
las <- seq(0.02, 2, length=3)
gas <- seq(0.02, 2, length=2)
nus <- seq(0.02, 2, length=2)
mus <- seq(0.02, 2, length=3)

allpar <- expand.grid(la=las, mu=mus, nu=nus, ga=gas)
allpar <- allpar[allpar$mu < 10 * allpar$la,]
allpar <- allpar[allpar$la < 10 * allpar$mu,]

## to avoid nested applies, expand allpar here (and make matrix to avoid vector
## problems in apply)
allpar <- as.matrix(allpar[rep(1:nrow(allpar), each=nRep), ])
allpar <- rbind(allpar, allpar, allpar, allpar)

## loop over paramters
noReturn <- mclapply((nrow(allpar)/2 + 1):nrow(allpar), mc.cores=20, FUN=function(i) {
	par <- allpar[i, ]
	x <- evolve(par[1], par[2], par[3], par[4], nGen=nGen)
	
	if(!is.null(nrow(x[[1]]))) {
		if(any(colSums(x[[1]], na.rm=TRUE) != 0)) {
			save(x, file=sprintf('constantRatesOut/simOut_constantRate_%s.RData', i))
		}
	}
})

