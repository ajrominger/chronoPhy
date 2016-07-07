setwd('~/Dropbox/Research/chronoPhy/poster')

library(ape)

## demonstration tree
tre <- rcoal(10)

pdf(file='fig_examplePhylo.pdf', width=3, height=2)
par(mar=rep(0, 4))
plot(tre, show.tip.label=FALSE, y.lim=c(1, 12))
segments(x0=max(branching.times(tre)) - sort(branching.times(tre))[-9], 
         y0=0, y1=10.5, col=hsv(alpha=0.5))
segments(x0=max(branching.times(tre)) - sort(branching.times(tre))[-9], 
         x1=0.1,
         y0=10.5, y1=12, col=hsv(alpha=0.5))
dev.off()


## gamma

pdf(file='fig_hiGammaPhylo.pdf', width=5, height=1.75)
par(mar=rep(0, 4))
tre1 <- read.tree(text='(a:2, (b:0.2, (c:0.1, d:0.1):0.1):1.8);')
plot(tre1, show.tip.label=FALSE, edge.width=3)
dev.off()

pdf(file='fig_loGammaPhylo.pdf', width=5, height=1.75)
par(mar=rep(0, 4))
tre2 <- read.tree(text='(a:2, (b:1.8, (c:1.6, d:1.6):0.2):0.2);')
plot(tre2, show.tip.label=FALSE, edge.width=3)
dev.off()

## chrono map