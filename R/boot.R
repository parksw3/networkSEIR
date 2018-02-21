correct.estfun <- function(gen, r) {
    wgen <- weighted.mean(gen, exp(r*gen))
    wR <- mean(exp(r*gen))
    
    return(c(gen=wgen, R=wR, r=r))
}

uncorrect.estfun <- function(gen, r) {
    mgen <- mean(gen)
    mR <- 1/mean(exp(-r*gen))
    
    return(c(gen=mgen, R=mR, r=r))
}

generation.bootstrap <- function(data, 
                                 level=0.95, 
                                 nsim=1000) {
    l <- (1-level)/2
    
    gen <- data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(data$day))
    
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
    
    bootlist <- replicate(nsim, sample(gen, length(gen), replace=TRUE), simplify=FALSE)
    
    boottest <- lapply(bootlist, function(x){
        r <- rnorm(1, mean=r, sd=sd.r)
        
        list(
            correct=correct.estfun(x, r),
            uncorrect=uncorrect.estfun(x, r)
        )
    })
    
    df <- boottest %>%
        lapply(bind_rows) %>%
        lapply(mutate, type=c("gen", "R", "r")) %>%
        bind_rows %>%
        gather(key, value, -type) %>%
        group_by(type, key) %>%
        summarize(
            lwr=quantile(value, 0.025),
            upr=quantile(value, 0.975)
        ) %>%
        arrange(key) %>%
        group_by %>%
        mutate(
            estimate=c(correct.estfun(gen, r)[c(1, 3, 2)], uncorrect.estfun(gen, r)[c(1, 3, 2)])
        )

    return(df)
}
