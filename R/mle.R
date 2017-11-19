gamma.R <- function(shape, scale, tau, r) (1+r*scale)^(shape)

censor.nll <- function(shape, scale,
                       tau, r) {
    R <- gamma.R(shape, scale, tau, r)
    
    -sum(dgamma(tau, shape=shape, scale=scale, log=TRUE) -r*tau + log(R))
}

censor.mle <- function(tau, r) {
    w <- exp(r*tau)
    w <- w/sum(w)
    
    mean <- weighted.mean(tau, w=w)
    var <- sum(w * (tau-mean)^2)
    
    shape <- mean^2/var
    
    scale <- mean/shape
    
    suppressWarnings(
        m <- bbmle::mle2(censor.nll
            , start=list(shape=shape, scale=scale)
            , data=list(tau=tau, r=r))
    )
    m
}

mle.R <- function(tau, r) {
    m <- censor.mle(tau, r)
    par <- bbmle::coef(m)
    unname(gamma.R(par[1], par[2], tau, r))
}
