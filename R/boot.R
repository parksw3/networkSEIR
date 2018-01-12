boot <- function(x, weight, level=0.95, nsim=1000) {
    n <- length(x)
    ss <- replicate(nsim,{
        s <- sample(1:n, n, replace=TRUE)
        weighted.mean(x[s], weight[s])
    })
    alpha <- (1-level)/2
    quantile(ss, c(alpha, 1-alpha))
}

boot_weight <- function(x, weight, level=0.95, nsim=1000) {
    n <- length(x)
    ss <- replicate(nsim, mean(sample(x, n, replace=TRUE, prob=weight)))
    alpha <- (1-level)/2
    quantile(ss, c(alpha, 1-alpha))
}

boot_half <- function(x, weight, level=0.95, nsim=1000) {
    n <- length(x)
    ss <- replicate(nsim,{
        s <- runif(n)<0.5
        weighted.mean(x[s], weight[s])
    })
    alpha <- (1-level)/2
    quantile(ss, c(alpha, 1-alpha), na.rm=TRUE)
}
