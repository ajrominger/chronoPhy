## script to run constant rate model over different par for 500 gen

setwd('~/Dropbox/Research/chronoPhy')
# source('chronoPhy_nonLinear3.R')
source('analysis_functions.R')
library(ape)
library(parallel)


## files to be analyzed
files <- list.files('expLaMuOut')

## only keep those that meet minimum completion requirements
keepem <- sapply(files, FUN=function(f) {
	load(paste('expLaMuOut', f, sep='/'))

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
	load(paste('expLaMuOut', f, sep='/'))
	island.summary(x)
}))


mean.stats.exp <- apply(island.stat[, c(1:3, 5:7), ], 2, function(x) apply(x, 1, mean, na.rm=TRUE))
se.stats.exp <- apply(island.stat[, 1:3, ], 2, function(x) apply(x, 1, sd, na.rm=TRUE))#/sqrt(sum(keepem))
