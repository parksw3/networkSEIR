source("ebola_param_seminr.R")
source("../R/seminr.R")
source("../R/generation.R")

N <- 10000

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    rr <- seminr.full(N, beta/(N-1), sigma, m.sigma, gamma, n.gamma, I0 = 10)
    if(nrow(rr$data) > 900) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="seminr_full_sim.rda")
