## assuming gamma... want non-parameteric eventually...
minuslogl.conditional <- function(log.mean, log.shape, data, tmax) {
    mean <- exp(log.mean)
    shape <- exp(log.shape)
    scale <- mean/shape
    
    -sum(
        dgamma(data$generation, shape=shape, scale=scale, log=TRUE) -
            pgamma(tmax-data$t_infected[match(data$infected_by, data$index)], shape=shape, scale=scale, log=TRUE),
        na.rm=TRUE
    )
}

minuslogl.full <- function(log.R, log.mean, log.shape, data, tmax) {
    R <- exp(log.R)
    mean <- exp(log.mean)
    shape <- exp(log.shape)
    scale <- mean/shape
    
    t_censor <- tmax-data$t_infected
    
    -(sum(
        log.R +
            dgamma(data$generation, shape=shape, scale=scale, log=TRUE),
        na.rm=TRUE
    ) + sum(
        -R * pgamma(t_censor, shape=shape, scale=scale),
        na.rm=TRUE
    ))
}

minuslogl.r <- function(log.lambda, log.r) {
    r <- exp(log.r)
    sum(
        log.lambda +
            r * data$t_infected
    ) + sum(
        - lambda * int_0^T 
    )
}

