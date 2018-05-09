library(bbmle)
source("../sim/full_param2.R")
source("../R/generation.R")
source("../R/empirical.R")
source("../R/boot.R")

load("../data/full_sim2.rda")

true.mean <- 1/sigma + 1/gamma

set.seed(101)
nQuant <- 10000
q <- (2*(1:nQuant)-1)/(2*nQuant)

lat <- qexp(q, rate=sigma)
inf <- qexp(q, rate=gamma)

ii <- sample(inf, prob=inf, replace=TRUE)

gen <- lat + runif(nQuant, min=0, max=ii)

## approximately true variance because I'm lazy...
true.var <- var(gen)

## assuming gamma unfortunately... want non-parameteric eventually...
minuslogl <- function(log.mean, log.var, data, tmax) {
    mean <- exp(log.mean)
    var <- exp(log.var)
    shape <- mean^2/var
    scale <- mean/shape

    -sum(
        dgamma(data$generation, shape=shape, scale=scale, log=TRUE) -
            pgamma(tmax-data$t_infected[match(data$infected_by, data$index)], shape=shape, scale=scale, log=TRUE),
        na.rm=TRUE
    )
}

rlist <- vector('list', 100)

set.seed(101)
for (i in 1:100) {
    cat(i)
    sim <- reslist[[i]]
    
    tmax <- max(sim$data$time)
    
    data <- generation.data(sim, tmax=tmax)
    
    data <- data[1:2000,]
    
    tmax <- max(data$t_infected)
    
    m <- mle2(minuslogl, 
        start=list(log.mean=log(true.mean), log.var=log(100)),
        method="Nelder-Mead",
        optimizer="optim",
        data=list(data=data, tmax=tmax)
    )
    
    cc <- coef(m)
    ci <- confint(m)
    
    res <- data.frame(
        param=c("mean", "var"),
        mean=exp(cc),
        lwr=exp(ci)[,1],
        upr=exp(ci)[,2],
        cover=c(
            exp(ci)[1,1] < true.mean && true.mean < exp(ci)[1,2],
            exp(ci)[2,1] < true.var && true.var < exp(ci)[2,2]
        )
    )
    
    print(res)
    
    rownames(res) <- NULL
    
    rlist[[i]] <- res
}

save("rlist", file="mle.rda")
