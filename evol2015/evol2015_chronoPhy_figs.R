setwd('~/Dropbox/Research/chronoPhy/evol2015')


library(maps)
library(mapdata)
library(ape)

## hawaii map
pdf(file='fig_hiMap.pdf', width=6, height=6)
map('worldHires', xlim=c(-160, -154), ylim=c(18, 23), res=0, fill=TRUE)
dev.off()

## trees
tre <- rphylo(50, 1, 0.8)

pdf(file='fig_bigTree.pdf', width=6, height=6)
plot(tre, show.tip.label=FALSE, edge.width=1.5)
dev.off()

tre.ka <- rphylo(10, 1, 0.8, fossils=TRUE)
tre.oa <- rphylo(8, 1, 0.5, fossils=TRUE)
tre.ma <- rphylo(6, 1, 0.5)
tre.ha <- rphylo(4, 1, 0.3)

pdf(file='fig_treeKa.pdf', width=3, height=3)
par(mar=rep(0, 4))
plot(tre.ka, x.lim=c(0, 6), show.tip.label=FALSE, edge.width=1.5)
dev.off()

pdf(file='fig_treeOa.pdf', width=3, height=3)
par(mar=rep(0, 4))
plot(tre.oa, x.lim=c(0, 6), y.lim=c(1, 20), show.tip.label=FALSE, edge.width=1.5)
dev.off()

pdf(file='fig_treeMa.pdf', width=3, height=3)
par(mar=rep(0, 4))
plot(tre.ma, x.lim=c(0, 6), y.lim=c(1, 20), show.tip.label=FALSE, edge.width=1.5)
dev.off()

pdf(file='fig_treeHa.pdf', width=3, height=3)
par(mar=rep(0, 4))
plot(tre.ha, x.lim=c(0, 6), y.lim=c(1, 20), show.tip.label=FALSE, edge.width=1.5)
dev.off()


## sampling distributions of extincton

x <- c(rnorm(140, 6, 1), abs(rnorm(100, 0, 0.1)))
x.den <- density(x, kernel='r', from=0)
x.den$x <- x.den$x - 1.5
x.den$y[x.den$x < 0] <- 0

pdf(file='fig_simExtSamp.pdf', width=4, height=4)
par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
plot(x.den, xlab='', xaxt='n', yaxt='n', main='', cex.lab=1.4, lwd=2)
abline(v=6, col='red', lwd=2)
mtext(expression(mu), side=1, line=2, cex=2)
axis(1, at=0)
dev.off()

x <- rgamma(100, 0.2, 5)
x.den <- density(x, kernel='r', from=0)

pdf(file='fig_obsExtSamp.pdf', width=4, height=4)
par(mar=c(3, 3, 0, 0)+0.1, mgp=c(1, 1, 0))
plot(x.den, xlab='', xaxt='n', yaxt='n', main='', cex.lab=1.4, lwd=2)
mtext(expression(mu), side=1, line=2, cex=2)
axis(1, at=0)
dev.off()


## small tree
tre <- rphylo(8, 1, 0.8)
pdf(file='fig_smallTre.pdf', width=4, height=4)
plot(tre, show.tip.label=FALSE, edge.width=2)
dev.off()



## SADs

pdf(file='fig_fisherSAD.pdf', width=2, height=2)
par(mar=c(2, 2, 0, 0), mgp=c(1, 1, 0))
plot(sort(rnbinom(300, mu=2, size=1) + 1, TRUE), log='y', xaxt='n', yaxt='n', xlab='Species rank', ylab='Abundance')
dev.off()

pdf(file='fig_lnormSAD.pdf', width=2, height=2)
par(mar=c(2, 2, 0, 0), mgp=c(1, 1, 0))
plot(sort(rpois(300, 2) + 1, TRUE), log='y', xaxt='n', yaxt='n', xlab='Species rank', ylab='Abundance')
dev.off()