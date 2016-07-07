## script to run constant rate model over different par for 500 gen

setwd('~/Dropbox/Research/chronoPhy')
# source('chronoPhy_nonLinear3.R')
source('analysis_functions.R')
library(ape)
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
# unique(allpar)

## files to be analyzed
files <- list.files('constantRatesOut')
parindex <- as.numeric(gsub('\\.RData', '', gsub('.*_', '', files)))

la.col <- hsv(0.7, 0.8, 0.7)
mu.col <- hsv(0.1, 0.4)
nu.col <- hsv(0.5, 0.5, 0.7)

if(!exists('posovershoot')) posovershoot <- FALSE

if(posovershoot) {
	## for possible overshoot use this:
	files <- files[parindex %in% which(allpar[, 1] == 2 & allpar[, 2] == 2 & 
	                                   allpar[, 3] == 2)]
	
	pdf(file='poster/fig_flatRatesHiNu.pdf', width=4, height=4)
	par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
	plot(1, type='n', ylim=c(0, 2.4),
	     axes=FALSE, frame.plot=TRUE, xlab='', ylab='')
	abline(h=c(1.85, 1.8, 2.3, 0.8), col=c(la.col, mu.col, nu.col, la.col), lwd=3, lty=c(1, 1, 1, 2))
	dev.off()
	
} else {
	files <- files[parindex %in% which(allpar[, 1] == 2 & allpar[, 2] == 1.01 & 
                                       allpar[, 3] == 0.02 & allpar[, 4] == 2)]
	
	pdf(file='poster/fig_flatRatesLoNu.pdf', width=4, height=4)
	par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
	plot(1, type='n', ylim=c(0, 2.2),
	     axes=FALSE, frame.plot=TRUE, xlab='', ylab='')
	abline(h=c(1.95, 1, 0.05, 0.8), col=c(la.col, mu.col, nu.col, la.col), lwd=3, lty=c(1, 1, 1, 2))
	dev.off()
}

## only keep those that meet minimum completion requirements
keepem <- sapply(files, FUN=function(f) {
	load(paste('constantRatesOut', f, sep='/'))

	endtime <- max(x[[3]]) > 4.5
	island.div <- sum(colSums(x[[1]]) > 0, na.rm=TRUE) > 0
	total.div <- sum(rowSums(x[[1]], na.rm=TRUE) > 0) > 8
	
	return(endtime & island.div & total.div)
})

sum(keepem)

files <- files[keepem]

## summarize stats across islands
## first position is diversity
## second is gamma
## remaining 5 are 0, 0.25, 0.5, 0.75, 1 quantiles of branching times
island.stat <- simplify2array(mclapply(files, mc.cores=3, FUN=function(f) {
	load(paste('constantRatesOut', f, sep='/'))
	island.summary(x)
}))

mean.stats.const <- apply(island.stat[, c(1:3, 5:7), ], 2, function(x) apply(x, 1, mean, na.rm=TRUE))
se.stats.const <- apply(island.stat[, 1:3, ], 2, function(x) apply(x, 1, sd, na.rm=TRUE))#/sqrt(sum(keepem))
mean.stats.const
se.stats.const