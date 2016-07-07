setwd('~/Dropbox/Research/chronoPhy')

pdf(width=4, height=4, file='fig_hypDiv.pdf')
par(mar=c(2, 2, 3, 0.2)+0.1, xpd=NA, mgp=c(1, 1, 0))
curve(dgamma(x, 2, 3), from=0, to=3, col='blue', lwd=3,
      axes=FALSE, frame.plot=TRUE, 
      xlab='Island age', ylab='Per capita rate',
      panel.first=segments(x0=seq(0.1, 2.7, length=5),
                           y0=par('usr')[3], y1=1*par('usr')[4], 
                           col='gray80'))
text(0.77, dgamma(0.57, 2, 3), labels=expression(lambda), 
     cex=2, col='blue')
     
text(seq(0.1, 2.7, length=5), rep(1.05*par('usr')[4], 5), 
     labels=paste('Island', 1:5), srt=45, adj=c(0, 1))
par(xpd=FALSE)

curve(exp(0.25*x) - 1, add=TRUE, col='red', lwd=3)
text(2.4, exp(0.25*2.6)-1, labels=expression(mu), 
     cex=2, col='red')

abline(h=0.25, col='blue', lty=2, lwd=3)
text(2.4, 0.32, labels=expression(gamma), 
     cex=2, col='blue')

abline(h=0.45, col='gray40', lwd=3)
text(2.4, 0.5, labels=expression(nu), 
     cex=2, col='gray40')

dev.off()