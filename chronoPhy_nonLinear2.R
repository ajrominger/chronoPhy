## *****************************
## Rate functions.

## In following functions, capital greeks are system-wide rates while
## lower case greeks are per capita rates. S refers to the species by
## island matrix. Per captia rates will be nSite long. functions
## return a vector either nSite long or (nSite^2 - nSite) long
## *****************************

library(ape)
library(spatstat)

source('~/Dropbox/Research/chronoPhy/evol2phylo.R')

## one way to get inhomogeneous birth-death
# fun <- AbscontDistribution(d=function(x) exp(-2 * (1 + 4*x - exp(-1/2*x))), low=0, low1=0)
# curve(fun@d(x), from=0, to=5)
# fun@q((1:6)/7)

## **********************************
## Linear rate functions
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

## non-linear rate functions (can be used for either La or Mu)
## in general these functions return a vector of sojour times, select event by
## picking the smallest sojourn time (i.e. the one that came first)

RhoExpR <- function(rho0, b, S, t0) {
	rho0 <- rho0 * colSums(S, na.rm=TRUE)
	
	mapply(these.t=t0, these.rho=rho0, FUN=function(these.t, these.rho) {
		# browser()
		if(these.rho==0 | these.t < 0) {
			out <- NA
		} else {
			ratefun <- function(x, y=0) {
				res <- these.rho * exp(-b*(x-these.t))
				res[x - these.t < 0] <- 0
				
				if(b < 0) res[x - these.t > 10] <- 0
				
				return(res)
			}
			
			## inhomogeneous poisson process on line segment 0:50 (50 because unlikely to go past that)
			# if(b < 0) browser()
			pp <- rpoisppOnLines(ratefun, 
			                     L=psp(x0=0, y0=0, x1=50, y1=0, window=owin(c(-1, 60), c(-1, 1)), 
			                                                                check=FALSE),
			                           # lmax=max(ratefun(c(0, 50)))
			                           lmax=NULL
			                           )
			if(pp$n==0) {
				out <- 60
			} else {
				out <- min(pp$x) - these.t
			}
		}
		
		return(out)
	})
}

# S1 <- matrix(0, nrow=10, ncol=5)
# for(i in 1:5) S1[1:(5-i+1), i] <- 1

# S2 <- S1
# S2[, 5] <- NA

# S3 <- matrix(0, nrow=10, ncol=5)
# S3[1, ] <- 1
# # S3[, 2:5] <- NA

# x <- replicate(500, mean(RhoExpR(1, 0, S3, 4:0), na.rm=TRUE))
# y <- replicate(500, mean(rexp(5, 1)))

# plot(density(x, na.rm=TRUE))
# lines(density(y), col='red')
## linear rate functions

NuR <- function(nu0, S, combo) {
	these.rates <- Nu(nu0, S, combo)
	out <- these.rates
	out[these.rates > 0] <- rexp(sum(these.rates > 0), these.rates[these.rates > 0])
	out[these.rates == 0] <- NA
	
	return(out)
}

GammaR <- function(ga0, S) {
	these.rates <- Gamma(ga0, S)
	out <- these.rates
	out[these.rates > 0] <- rexp(sum(these.rates > 0), these.rates[these.rates > 0])
	out[these.rates == 0] <- NA
	
	return(out)
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

## sample from event space (... is for arguments to LaFun)
pickEvent <- function(LaFun, mu0, bmu, nu0, ga0, S, time, island.time, combo, nSite, ...) {
	t0 <- time - island.time
	
	## for islands that have not yet emerged we want rate to be 0, this is achieved by having elements of
	## S be NA if the island does not exist yet
	theseTimes <- list(L = LaFun(S=S, t0=t0, ...), M = RhoExpR(mu0, bmu, S, t0), G = GammaR(ga0, S), N = NuR(nu0, S, combo))
	event.index <- which.min(unlist(theseTimes))
	thisTime <- unlist(theseTimes)[event.index]
	
	if(length(event.index) == 0) browser()
	
	if(event.index <= nSite) { # lambda
		thisType <- 1
	} else if(event.index <= nSite*2) { # mu
		thisType <- 2
		event.index <- event.index - nSite
	} else if(event.index <= nSite*3) { # gamma
		thisType <- 4
		event.index <- event.index - nSite*2
	} else { # nu
		thisType <- 3
		event.index <- event.index - nSite*3
	}
	
	return(list(time=thisTime, 
	            type=thisType,
	            event.index=event.index))
}

## update functions for each type of event
event <- function(time, type, event.index) {
	if(type == 3) { # immigration
		thisImm <- fromTo(event.index, nsite = ncol(S))
		spp <- sample(which(S[, thisImm[1]] > 0 & S[, thisImm[2]] == 0), 1)
		
		S[spp, thisImm[2]] <<- 1
		cat('immigration \n')
	} else { # not immigration
		## site given by pickEvent
		site <- event.index
		
		## choose parent (or sp to go extinct) at random from w/n site
		parent <- sample(which(S[, site] > 0), 1)

		if (type %in% c(1, 4)) { # speciation
			cat('speciation')
			## new species is just next row in S
			newSp <- max(which(p > 0)) + 1

			## update S, p, x (time; time is being updated at each iteration so only need
			## to record it for speciation events)
			S[newSp, site] <<- 1
			p[newSp] <<- parent
			x[newSp] <<- time

			## if anagenic speciation also remove ancestor at site
			if(type == 4) {
				cat(' (anagenic)')
				S[parent, site] <<- 0
			}
			cat('\n')
		} else { # extinction
			cat('extinction')
			S[parent, site] <<- 0
		}
	}
	
	return(NULL)
}



## *****************************
## master function that evolves the system
## *****************************

evolve <- function(LaFun, mu0, bmu, nu0, ga0, nGen=100, ...) { # ... for LaFun
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
		print(paste(i, ': t =', times[i]))
		
		## update S based on what islands are emerged at begining of time period
		## put 0 (i.e. availible) where there were NA
		S[, siteAge <= times[i] & is.na(S[1, ])] <- 0
		
		## if all extinct then don't update anymore
		if(all(colSums(S, na.rm=TRUE) == 0)) {
			cat('all extinct \n')
			break
		}
		
		## pick event type and calcualte sojurn time
		thisEvent <- pickEvent(LaFun=LaFun, mu0=mu0, bmu=bmu, nu0=nu0, ga0=ga0, S=S, 
		                       time=times[i], island.time=siteAge, combo=combo, nSite=nSite, ...)
		
		## update time by sojurn
		times[i + 1] <- times[i] + thisEvent[[1]]
		thisEvent[[1]] <- times[i + 1]
		if(thisEvent[[1]] > 5) {
			cat('exceeded time \n')
			break
		}
			
		## update system with realization of the event
		event(thisEvent[[1]], thisEvent[[2]], thisEvent[[3]])
		
		cat('\n')
	}
	
	nspp <- max(which(p > 0))
	p[1] <- 0
	return(list(S[1:nspp, ], p[1:nspp], x[1:nspp]))
}


# x.phy <- evol2phylo(x, max(x[[3]])+1)
# plot(x.phy)

# abline(v=x[[3]][-1] - x[[3]][2], col=rgb(1, 0, 0, alpha=0.4))

# xmrk <- pretty(c(0, max(x[[3]][-1])))
# axis(1, at=xmrk, label=rev(xmrk))

# y <- evolve(la=0, mu=0, nu=1, ga=0, nGen=10)
