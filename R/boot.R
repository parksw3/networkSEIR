estfun <- function(gen, r) {
    wgen <- weighted.mean(gen, exp(r*gen))
    wR <- mean(exp(r*gen))
    
    return(c(gen=wgen, R=wR, r=r))
}

generation.bootstrap <- function(data, 
                                 level=0.95, 
                                 nsim=1000,
                                 true.r) {
    l <- (1-level)/2
    
    gen <- data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(data$day))
    
    if (missing(true.r)) {
        Sigma <- matrix(0, 3, 3)
        diag(Sigma) <- 0.0001
        
        start <- MASS::mvrnorm(100, mu=c(log.i0=log(0.1), r=0.05, log.size=log(30)), Sigma = Sigma)
        
        fitlist <- apply(start, 1, function(x) {
            fit <- suppressWarnings(
                bbmle::mle2(
                    case~dnbinom(mu=exp(log.i0)*exp(r*day), size=exp(log.size))
                    , start=as.list(x)
                    , data=data.frame(
                        case=case
                        , day=sort(unique(data$day))
                    )
                    , method="BFGS"
                )
            )
        })
        
        mle <- fitlist[[which.max(sapply(fitlist, logLik))]]
        
        r <- unname(coef(mle)[2])
        
        sd.r <- stdEr(mle)[2]
    } else {
        r <- true.r
    }
    
    bootlist <- replicate(nsim, sample(gen, length(gen), replace=TRUE), simplify=FALSE)
    
    bootest <- sapply(bootlist, function(x, true.r){
        if (missing(true.r)) {
            r <- rnorm(1, mean=r, sd=sd.r)
        } else {
            r <- true.r
        }
        estfun(x, r)
    }, true.r=true.r)
    
    res <- apply(bootest, 1, quantile, c(l, 1-l))
    estimate <- estfun(gen, r)
    
    df <- data.frame(
        type=c("gen", "R", "r"),
        estimate=estimate,
        lwr=res[1,],
        upr=res[2,]
    )
    
    rownames(df) <- NULL
    
    return(df)
}
