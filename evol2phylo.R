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
        
        
        # if(sum(edge[rows, 2] > Ntip) == 1) {
        	# edge.length[rows[edge[rows, 2] > Ntip]] <- s[child + 1] - s[child]
        # } else if(sum(edge[rows, 2] > Ntip) > 1) {
        	# new.kid1 <- Ntip - ceiling(which(edge[, 1] == edge[rows[1], 2]) / 2)[1] + 1
        	# new.kid2 <- Ntip - ceiling(which(edge[, 1] == edge[rows[2], 2]) / 2)[1] + 1
        	
        	# edge.length[rows[1]] <- s[new.kid1] - s[child]
        	# edge.length[rows[2]] <- s[new.kid2] - s[child]
        # }
         
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


# x <- evolve(1, 0, 0, 0, 10)
# x.phy <- evol2phylo(x)
# par(mar=c(5, 0, 0, 0))
# plot(x.phy)
# abline(v=c(x[[3]]) - x[[3]][2], col=hsv(alpha=0.5))
# axis(1, at=x[[3]][-1]-x[[3]][2], labels=round(x[[3]][-1], 4), las=2)


