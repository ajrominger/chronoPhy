## script to run constant rate model over different par for 500 gen

setwd('~/Dropbox/Research/chronoPhy')
source('chronoPhy.R')
source('analysis_functions.R')
library(ape)
library(parallel)

## set up simulation par
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

## files to be analyzed
files <- list.files('constantRatesOut')
par.index <- as.numeric(gsub('\\..*', '', gsub('.*_', '', files)))

## only keep those that meet minimum completion requirements
keepem <- sapply(files, FUN=function(f) {
	load(paste('constantRatesOut', f, sep='/'))

	endtime <- max(x[[3]]) > 4
	# island.div <- sum(colSums(x[[1]]) > 0, na.rm=TRUE) > 3
	total.div <- sum(rowSums(x[[1]], na.rm=TRUE) > 0) > 8
	
	return(endtime & island.div & total.div)
})

sum(keepem)

files <- files[keepem]
par.index <- par.index[keepem]
allpar <- allpar[par.index, ]

## summarize stats across islands
## first position is diversity
## second is gamma
## remaining 5 are 0, 0.25, 0.5, 0.75, 1 quantiles of branching times
island.stat <- simplify2array(mclapply(files, mc.cores=3, FUN=function(f) {
	load(paste('constantRatesOut', f, sep='/'))
	island.summary(x)
}))

## break up simulation runs by unique parameter combos and summarize stats accordingly
allpar <- as.data.frame(allpar)

allpar$growth <- log(allpar$la / allpar$mu)
allpar$growthclass <- cut(allpar$growth, 3)
allpar$immclass <- cut(allpar$nu, 2)


parunique <- unique(allpar[, c('growth', 'nu', 'growthclass', 'immclass')])
parcombo <- do.call(paste, allpar[, c('growthclass', 'immclass')])
table(parcombo)


chrono.div <- sapply(split(as.data.frame(t(island.stat[, 1, ])), parcombo), 
                     function(x) apply(x, 2, mean, na.rm=TRUE))
chrono.gam <- sapply(split(as.data.frame(t(island.stat[, 2, ])), parcombo), 
                     function(x) apply(x, 2, mean, na.rm=TRUE))
nrow(parunique)
ncol(chrono.div)


ylim <- quantile(as.vector(chrono.div), prob=c(0, 0.9))

palette(rainbow(ncol(chrono.div), end=0.8))

par(mfcol=c(3, ceiling(ncol(chrono.div)/3)), mar=rep(0, 4), oma=rep(1, 4))

for(i in 1:ncol(chrono.div)) {
	plot(chrono.div[, i], type='l', col=i, axes=FALSE)
	legend('topright', legend=i, cex=0.8, bty='n')
	box()
}

posOver <- scan('constantRatePossibleOvershoot.txt')
plot(density(log(parunique[, 1] / parunique[, 2])))
lines(density(log(parunique[posOver, 1] / parunique[posOver, 2])))


plot(density(parunique[, 3]))
lines(density(parunique[posOver, 3]))
