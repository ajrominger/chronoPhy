source('~/R_functions/logAxis.R')

## function to turn evolve output into ape phylo object

evol2phylo <- function(x, maxs) {
    p <- x[[2]]
    s <- x[[3]]
    if(missing(maxs)) maxs <- max(s)
    Ntip <- length(p)
    Nnode <- Ntip - 1
    tip <- 1:Ntip
    edge <- matrix(0, nrow=2*Ntip - 2, ncol=2)
    edge.length <- numeric(2*Ntip - 2)

    rows <- 1:2
    
    for(i in 1:(Ntip - 1)) {
        ## start from the end of the parent vector and work back
        child <- Ntip - i + 1
        
        edge[rows, 1] <- Ntip + i
        edge[rows, 2] <- c(tip[child], p[child])
        
        edge.length[rows[edge[rows, 2] <= Ntip]] <- maxs - s[child]
        
        if(edge[rows[1], 2] > Ntip) {
        	new.kid1 <- Ntip - ceiling(which(edge[, 1] == edge[rows[1], 2]) / 2)[1] + 1
        	edge.length[rows[1]] <- s[new.kid1] - s[child]
        }
        
        if(edge[rows[2], 2] > Ntip) {
        	new.kid2 <- Ntip - ceiling(which(edge[, 1] == edge[rows[2], 2]) / 2)[1] + 1
        	edge.length[rows[2]] <- s[new.kid2] - s[child]
        }
         
        tip[tip %in% c(child, p[child])] <- Ntip + i
        p[p %in% c(child, p[child])] <- Ntip + i
        
        rows <- rows + 2
    }
    
    
    
    out <- list(edge=edge, tip.label=1:Ntip, Nnode=Nnode, edge.length=edge.length)
    # browser()
    class(out) <- 'phylo'
    out <- read.tree(text=write.tree(out))

    return(out)
}


island.summary <- function(x) {
	phy <- evol2phylo(x)
	nIsland <- rowSums(x[[1]])
	
	out <- array(NA, dim=c(ncol(x[[1]]), 8))
	
	for(i in 1:ncol(x[[1]])) {
		these.extinct <- x[[1]][, i] == 0
		
		endemic <- sum(x[[1]][, i] == 1 & nIsland == 1)
		div <- sum(!these.extinct)
		
		if(div > 2) {
			this.phy <- drop.tip(phy, which(these.extinct))
			gam <- gammaStat(this.phy)
			bt <- quantile(branching.times(this.phy), prob=c(0, 0.25, 0.5, 0.75, 1))
			
			out[i, ] <- c(div, endemic, gam, bt)
		} else {
			out[i, ] <- c(div, endemic, rep(NA, 6))
		}
	}
	
	return(out)
}


plot.island.summary <- function(means, ses, ages=c(0.5, 1.5, 2.2, 3.8, 5), ylims, cols, ...) {
	# par(mfrow=c(1, 4), oma=c(3, 0, 0, 0)+0.1, mar=c(0, 3, 0, 0))
	
	if(missing(ylims)) ylims <- vector('list', 5)
	if(missing(cols)) cols <- rep('black', 5)
	poly.col <- rgb(t(col2rgb(cols)), maxColorValue=255, alpha=255/4)
	
	for(i in 1:4) {
		if(i < 4) {
			yse <- c(means[, i] + ses[, i], rev(means[, i] - ses[, i]))
			if(i < 3) yse[yse < 1] <- 1
		} else {
			yse <- c(means[, i], rev(means[, i+2]))
			i <- i + 1
		}
		
		if(is.null(ylims[[i]])) ylims[[i]] <- range(yse)
		
		plot(ages, means[, i], log=ifelse(i < 3, 'y', ''), type='b', 
		     yaxt=ifelse(i < 3, 'n', 's'), 
		     col=cols[ifelse(i < 4, i, i-1)],
		     panel.first = {
		     	polygon(x=c(ages, rev(ages)), y=yse, col=poly.col[ifelse(i < 4, i, i-1)], border=NA)
		     }, ylim=ylims[[i]], 
		     ylab=switch(i, 
		                 '1'='Total richness',
		                 '2'='Endemics',
		                 '3'='Gamma statistic',
		                 '4'='Branching',
		                 '5'='Branching'), cex.lab=1.5, ...)
		
		if(i < 3 & !('yaxt' %in% names(list(...)))) logAxis(2)
	}
}

# plot.island.summary(mean.stats.const, se.stats.const, cols=c('red', 'blue', 'green', 'black'))

