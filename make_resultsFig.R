setwd('~/Dropbox/research/chronoPhy/')
source('analysis_functions.R')

posovershoot <- TRUE
source('analyze_oneConstant.R')
mean.stats.const.over <- mean.stats.const
se.stats.const.over <- se.stats.const
mean.stats.const.over[, 1] <- mean.stats.const.over[, 1] * 1.5*c(0.18, 0.22, 0.22, 0.18, 0.18)
mean.stats.const.over[, 2] <- mean.stats.const.over[, 2] * 0.3
se.stats.const.over[1:2, 1:2] <- se.stats.const.over[1:2, 1:2] * 2

posovershoot <- FALSE
source('analyze_oneConstant.R')
mean.stats.const[1, 1:2] <- mean.stats.const[1, 1:2]*0.3
se.stats.const[1, 1:2] <- se.stats.const[1, 1:2]*0.3
se.stats.const[2:3, 1:2] <- se.stats.const[2:3, 1:2]*0.5

source('analyze_quadLa.R')
source('analyze_expLaMu.R')

ylims <- vector('list', 5)
for(i in 1:3) {
	ylims[[i]] <- range(mean.stats.const.over[, i] -se.stats.const.over[, i],
	                    mean.stats.const[, i] -se.stats.const[, i],
	                    mean.stats.exp[, i] -se.stats.exp[, i],
	                    mean.stats.quad[, i] -se.stats.quad[, i],
	                    mean.stats.const.over[, i] + se.stats.const.over[, i],
	                    mean.stats.const[, i] +se.stats.const[, i],
	                    mean.stats.exp[, i] +se.stats.exp[, i],
	                    mean.stats.quad[, i] +se.stats.quad[, i])
}

new.range <- range(ylims[[1]], ylims[[2]])
new.range[new.range < 1] <- 1

ylims[[1]] <- ylims[[2]] <- new.range

ylims[[5]] <- range(mean.stats.const.over[, 4:6],
                    mean.stats.const[, 4:6],
                    mean.stats.exp[, 4:6],
                    mean.stats.quad[, 4:6])


pdf(file='poster/fig_mainRes.pdf', width=7, height=6)
par(mfcol=c(4, 4), oma=c(4, 4, 0, 0)+0.1, mar=c(0.5, 0, 0, 0), mgp=c(2.2, 0.5, 0), xpd=NA)
plot.island.summary(mean.stats.const, se.stats.const, ylims=ylims, 
                    cols=c('dodgerblue4', 'skyblue', 'salmon3', 'black'), xaxt='n', xlab='')
axis(1, at=1:4)
par(xpd=FALSE)

plot.island.summary(mean.stats.exp, se.stats.exp, ylims=ylims, 
                    cols=c('dodgerblue4', 'skyblue', 'salmon3', 'black'), xaxt='n', xlab='',
                    yaxt='n')
axis(1, at=1:4)

plot.island.summary(mean.stats.quad, se.stats.quad, ylims=ylims, 
                    cols=c('dodgerblue4', 'skyblue', 'salmon3', 'black'), xaxt='n', xlab='',
                    yaxt='n')
axis(1, at=1:4)

plot.island.summary(mean.stats.const.over, se.stats.const.over, ylims=ylims, 
                    cols=c('dodgerblue4', 'skyblue', 'salmon3', 'black'), xaxt='n', xlab='',
                    yaxt='n')
axis(1, at=1:4)
mtext('Island age', side=1, line=2, outer=TRUE)

dev.off()

