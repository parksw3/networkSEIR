estfun <- function(data) {
    gen <- data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(data$day))
    
    fit <- suppressWarnings(
        bbmle::mle2(case~dpois(lambda=i0*exp(r*day)),
                    start=list(i0=1, r=0.05), 
                    data=data.frame(
                        case=case,
                        day=sort(unique(data$day)),
                        method="BFGS"
                        ))
    )
    
    r <- unname(coef(fit)[2])
    wgen <- weighted.mean(gen, exp(r*gen))
    wR <- mean(exp(r*gen))
    c(gen=wgen, R=wR, r=r)
}

cifun <- function(data) {
    gen <- data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(data$day))
    
    fit <- suppressWarnings(
        bbmle::mle2(case~dpois(lambda=i0*exp(r*day)),
                    start=list(i0=1, r=0.05), 
                    data=data.frame(
                        case=case,
                        day=sort(unique(data$day)),
                        method="BFGS"
                    ))
    )
    
    r <- unname(coef(fit)[2])
    wgen <- weighted.mean(gen, exp(r*gen))
    wR <- mean(exp(r*gen))
    c(gen=wgen, R=wR, r=r)
}

bootfun <- function(data, method=c("resample", "half", "jackknife"), i) {
    method <- match.arg(method)
    
    if (method=="resample") {
        index <- sample(1:nrow(data), nrow(data), replace=TRUE)
    } else if (method=="half") {
        index <- (1:nrow(data))[runif(nrow(data)) < 0.5]
    } else if (method=="jackknife") {
        index <- (1:nrow(data))[-i]
    }
    
    filter.data <- data[index,]
    
    estfun(filter.data)
}

generation.bootstrap <- function(data, 
                                 level=0.95, 
                                 nsim=1000, 
                                 method=c("resample", "half"),
                                 bca=FALSE) {
    l <- (1-level)/2

    estimate <- estfun(data)
    
    bsample <- replicate(nsim,bootfun(data, method))
    
    res <- apply(bsample, 1, quantile, c(l, 1-l), na.rm=TRUE)
    
    rownames(res) <- NULL
    
    data.frame(
        type=c("gen", "R", "r"),
        estimate=estimate,
        lwr=res[1,],
        upr=res[2,]
    )
}
