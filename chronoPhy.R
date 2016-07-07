## *****************************
## Rate functions.

## In following functions, capital greeks are system-wide rates while
## lower case greeks are per capita rates. S refers to the species by
## island matrix. Per captia rates will be nSite long. functions
## return a vector either nSite long or (nSite^2 - nSite) long
## *****************************

library(ape)
library(distr)

source('~/Dropbox/Research/chronoPhy/evol2phylo.R')

## one way to get inhomogeneous birth-death
# fun <- AbscontDistribution(d=function(x) exp(-2 * (1 + 4*x - exp(-1/2*x))), low=0, low1=0)
# curve(fun@d(x), from=0, to=5)
# fun@q((1:6)/7)

## **********************************
## Rate functions:
## 
## **********************************

## cladogenic speciation
Lambda <- function(la, S) {
	la * colSums(S, na.rm=TRUE)
}

## island extinction
Mu <- function(mu, S) {
	mu * colSums(S, na.rm=TRUE)
}

## unique immigration
Nu <- function(nu, S, combo) {
	nu * nUnique(S, combo)
}

## anagenic speciation leading to increased S system-wide but not on an island
Gamma <- function(ga, S) {
	ga * nShared(S)
}

## helper function to find number of spp unique to i wrt j. if j is NA then rate is 0
nUnique <- function(S, combo) {
	apply(combo, 1, function(x) sum(S[, x[1]] == 1 & S[, x[2]] == 0 & !is.na(S[, x[2]]), na.rm=TRUE))
}

## helper function to find number of spp shared between i and any other j
nShared <- function(S) {
	## finds pairs of islands that have at least 1 sp shared
	couldAna <- rowSums(S, na.rm=TRUE) > 1
	
	## sums up number shared
	if(sum(couldAna) == 1) {
		out <- S[couldAna, ]
		out[is.na(out)] <- 0
		return(out)
	} else {
		out <- colSums(S[rowSums(S, na.rm=TRUE) > 1, ])
		out[is.na(out)] <- 0
		return(out)
	}
}


## *****************************
## Functions to check for which type of event occured and update
## system based on event type.

## S is species by island matrix, x is vector of times, p is vecotor
## of parent nodes
## *****************************

## helper function when determining where immigrants come from and go
fromTo <- function(n, nsite) {
	## n is a number in 1:(nsite^2) corresponding to the cell of an nsite x nsite
	## matrix so these two lines decide which origin site n refers to (i.e. the row)
	from <- n%%nsite
	from[from == 0] <- nsite
	
	## similarely ceiling(n/nsite) determines what destination site n referst to 
	## (i.e. the column)
	c(from = from, to = ceiling(n/nsite))
}

## sample from event space
pickEvent <- function(la, mu, nu, ga, S, x, p, combo) {
	# browser()
	## for islands that have not yet emerged we want rate to be 0, this is achieved by having elements of
	## S be NA if the island does not exist yet
	theseRates <- list(L = Lambda(la, S), M = Mu(mu, S), N = Nu(nu, S, combo), G = Gamma(ga, S))
	
	## everything extinct if rates are 0 (shouldn't ever come up)
	# if(all(unlist(theseRates) == 0)) {
		# cat('all extinct \n')
		# return(list(rho=0, time=0, type=5))
	# }
	
	## here the rates should be used like at begining of script (in example with distr)
	thisTime <- rexp(1, sum(unlist(theseRates)))
	
	## for the sampling of event types should the `rates' be integrals of rate functions over
	## sojourn time or should they be instantaneous rates at event time?
	thisType <- sample(4, 1, prob = sapply(theseRates, sum))
	
	return(list(rho=theseRates[[thisType]], 
	            time=thisTime, 
	            type=thisType))
}

## update functions for each type of event
event <- function(rho, time, type) {
	## if immigration need to modify rho s.t. refers to single site
	## rates
	if (type == 3) { # immigration
		thisImm <- fromTo(sample(length(rho), 1, prob = rho), nsite = ncol(S))
		spp <- sample(which(S[, thisImm[1]] > 0 & S[, thisImm[2]] == 0), 1)
		S[spp, thisImm[2]] <<- 1
		# cat(thisImm[1], 'to ', thisImm[2], '\n')
	} else { # not immigration
		## site choosen at random, weighted by rates
		site <- try(sample(ncol(S), 1, prob = rho))
		if(inherits(site, 'try-error')) browser()
		# cat(site, '\n')
		
		## choose parent (or sp to go extinct) at random from w/n site
		parent <- sample(which(S[, site] > 0), 1)

		if (type %in% c(1, 4)) { # speciation
			## new species is just next row in S
			newSp <- max(which(p > 0)) + 1

			## update S, p, x (time; time is being updated at each iteration so only need
			## to record it for speciation events)
			S[newSp, site] <<- 1
			p[newSp] <<- parent
			x[newSp] <<- time

			## if anagenic speciation also remove ancestor at site
			if (type == 4) {
				S[parent, site] <<- 0
			}
		} else { # extinction
			S[parent, site] <<- 0
		}
	}
	
	return(NULL)
}



## *****************************
## master function that evolves the system
## *****************************

evolve <- function(la, mu, nu, ga, nGen=100) {
	## make event see objects in internal environment
	environment(event) <- environment()
	
	## make objects acted on by event
	nMax <- 10^5
	nSite <- 5
	siteAge <- 5 - c(5, 2.5, 1.5, 1.2, 0.5)
	combo <- expand.grid(1:nSite, 1:nSite)
	
	## the key objects that get updated by `event'
	S <- matrix(0, nrow = nMax, ncol = nSite)
	p <- numeric(nMax)
	x <- numeric(nMax)
	
	## initialize with 1 sp in first site at time 0 and sites not yet formed as NA
	S[1, 1] <- 1
	p[1] <- 1
	S[, siteAge > max(p)] <- NA
	
	## keep track of total system time (first time is 0 = initial state)
	times <- numeric(nGen + 1)
	
	## can't be parallelized in a strait-forward way because step i+1 depends on i
	for(i in 1:nGen) {
		# cat('i = ', i, ': ')
		## update S based on what islands are emerged at begining of time period
		## put 0 (i.e. availible) where there were NA
		S[, siteAge <= times[i] & is.na(S[1, ])] <- 0
		
		# if(i > 1 & any((siteAge <= times[i]) != (siteAge <= times[i-1]))) {
			# cat('**island', which((siteAge <= times[i]) != (siteAge <= times[i-1])), 'formed \n')
		# }
		
		## if all extinct then don't update anymore
		if(all(colSums(S, na.rm=TRUE) == 0)) {
			# cat('all extinct \n')
			break
		}
		
		## pick event type and calcualte sojurn time
		thisEvent <- pickEvent(la, mu, nu, ga, S, x, p, combo)
		
		## update time by sojurn
		times[i + 1] <- times[i] + thisEvent[[2]]
		thisEvent[[2]] <- times[i + 1]
		if(thisEvent[[2]] > 5) {
			# cat('exceeded time \n')
			break
		}
		
		# cat(switch(thisEvent[[3]], 
		           # '1' = 'cladogenic speciation on island ',
		           # '2' = 'extinction on island ',
		           # '3' = 'immigration from island ',
		           # '4' = 'anagenic speciatoin on island '))
		
		## update system with realization of the event
		event(thisEvent[[1]], thisEvent[[2]], thisEvent[[3]])
		
		# print(S[p > 0, ])
	}
	
	nspp <- max(which(p > 0))
	p[1] <- 0
	return(list(S[1:nspp, ], p[1:nspp], x[1:nspp]))
}

x <- evolve(la=0.4, mu=0, nu=0, ga=0, nGen=10)
x

# x.phy <- evol2phylo(x, max(x[[3]])+1)
# plot(x.phy)

# abline(v=x[[3]][-1] - x[[3]][2], col=rgb(1, 0, 0, alpha=0.4))

# xmrk <- pretty(c(0, max(x[[3]][-1])))
# axis(1, at=xmrk, label=rev(xmrk))

# y <- evolve(la=0, mu=0, nu=1, ga=0, nGen=10)
