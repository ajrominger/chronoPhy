gamma.from.btime <- function(bt) {
    N <- length(bt)
    bt <- bt[-1] # remove the arbitrary 0
    g <- rev(c(bt[1], diff(bt)))
    ST <- sum((2:N) * g)
    stat <- sum(cumsum((2:(N - 1)) * g[-(N - 1)]))/(N - 2)
    m <- ST/2
    s <- ST * sqrt(1/(12 * (N - 2)))
    (stat - m)/s
}