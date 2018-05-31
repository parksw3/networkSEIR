source("ebola_param.R")
source("kernel.R")
source("../R/generation.R")

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    rr <- kernel.full(N, beta/gamma, function(x) rgamma(x, shape=true.shape, scale=true.mean/true.shape), I0 = 10, imax=1000)
    if(nrow(rr$data) > 900) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="kernel_sim.rda")
