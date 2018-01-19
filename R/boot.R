bootfun <- function(data, method=c("resample", "half")) {
    method <- match.arg(method)
    
    if (method=="resample") {
        index <- sample(1:nrow(data), nrow(data), replace=TRUE)
    } else {
        index <- (1:nrow(data))[runif(nrow(data)) < 0.5]
    }
    
    filter.data <- data[index,]
    
    gen <- filter.data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(filter.data$week))
    
    fit <- lm(log(case)~sort(unique(filter.data$week)))
    r <- coef(fit)[2]/7
    
    c(gen=weighted.mean(gen, exp(r*gen)), R=mean(exp(r*gen)), r=r)
}

generation.bootstrap <- function(data, level=0.95, nsim=1000, method=c("resample", "half")) {
    l <- (1-level)/2

    gen <- data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(data$week))
    
    fit <- lm(log(case)~sort(unique(data$week)))
    r <- coef(fit)[2]/7
    
    res <- apply(replicate(nsim,bootfun(data, method)), 1, quantile, c(l, 1-l), na.rm=TRUE)
    rownames(res) <- NULL
    
    data.frame(
        type=c("gen", "R", "r"),
        estimate=c(weighted.mean(gen, exp(r*gen)), mean(exp(r*gen)), r),
        lwr=res[1,],
        upr=res[2,]
    )
}
