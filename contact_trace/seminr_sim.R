source("ebola_param_seminr.R")
source("../R/seminr.R")
source("../R/generation.R")

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    rr <- seminr.full(N, beta/N, sigma, m.sigma, gamma, n.gamma, I0 = 10, imax=1000)
    if(nrow(rr$data) > 900) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="seminr_sim.rda")
