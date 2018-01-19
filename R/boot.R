bootfun <- function(data) {
    index <- (1:nrow(data))[runif(nrow(data)) < 0.5]
    filter.data <- data[index,]
    
    gen <- filter.data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(filter.data$week))
    
    fit <- lm(log(case)~unique(filter.data$week))
    r <- coef(fit)[2]/7
    
    c(gen=weighted.mean(gen, exp(r*gen)), R=mean(exp(r*gen)), r=r)
}

generation.bootstrap <- function(data, level=0.95, nsim=1000) {
    l <- (1-level)/2

    gen <- data$gen
    gen <- gen[!is.na(gen)]
    
    case <- as.vector(table(data$week))
    
    fit <- lm(log(case)~unique(data$week))
    r <- coef(fit)[2]/7
    
    res <- apply(replicate(nsim,bootfun(data)), 1, quantile, c(l, 1-l), na.rm=TRUE)
    
    data.frame(
        type=c("gen", "R", "r"),
        estimate=c(weighted.mean(gen, exp(r*gen)), mean(exp(r*gen)), r),
        lwr=res[1,],
        upr=res[2,]
    )
}
