library(bbmle)
source("../sim/gamma_param.R")
source("../R/generation.R")
source("../R/empirical.R")

load("../data/gamma_sim.rda")

true.mean <- 1/sigma + 1/gamma
true.shape <- 2

true.R <- beta/gamma

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

rlist <- vector('list', 100)

set.seed(101)
for (i in 1:100) {
    cat(i)
    sim <- reslist[[i]]
    
    tmax <- max(sim$data$time)
    
    data <- generation.data(sim, tmax=tmax)[1:1000,]
    
    tmax <- max(data$t_infected)
    
    m1 <- mle2(minuslogl.conditional, 
        start=list(log.mean=log(true.mean), log.shape=log(true.shape)),
        method="Nelder-Mead",
        optimizer="optim",
        data=list(data=data, tmax=tmax),
        control=list(maxit=10000)
    )
    
    pp1 <- profile(m1)
    
    cc1 <- coef(m1)
    ci1 <- confint(pp1)
    
    res1 <- data.frame(
        method="conditional",
        param=c("R", "mean", "shape"),
        mean=c(NA, exp(cc1)),
        upr=c(NA, exp(ci1)[,2]),
        lwr=c(NA, exp(ci1)[,1]),
        cover=c(
            NA,
            exp(ci1)[1,1] < true.mean && true.mean < exp(ci1)[1,2],
            exp(ci1)[2,1] < true.shape && true.shape < exp(ci1)[2,2]
        )
    )
    
    print(res1)
    
    m2 <- mle2(minuslogl.full, 
               start=list(log.R=log(true.R), log.mean=log(true.mean), log.shape=log(true.shape)),
               method="Nelder-Mead",
               optimizer="optim",
               data=list(data=data, tmax=tmax),
               control=list(maxit=10000)
    )
    
    pp2 <- profile(m2)
    
    cc2 <- coef(m2)
    ci2 <- confint(pp2)
    
    res2 <- data.frame(
        method="full",
        param=c("R", "mean", "shape"),
        mean=exp(cc2),
        upr=exp(ci2)[,2],
        lwr=exp(ci2)[,1],
        cover=c(
            exp(ci2)[,1] < c(true.R, true.mean, true.shape) &
                exp(ci2)[,2] > c(true.R, true.mean, true.shape)
        )
    )
    
    print(res2)
    
    rownames(res1) <- rownames(res2) <- NULL
    
    rlist[[i]] <- do.call("rbind", list(res1, res2))
}

save("rlist", file="gamma_mle.rda")
