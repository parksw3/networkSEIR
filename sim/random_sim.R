library(igraph)
source("cmp_param.R")
source("../R/seir.R")

load("../data/random_network.rda")

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    
    rr <- seir(g, beta, sigma, gamma, I0=10, imax=500)
    if(nrow(rr$data) > 400) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="random_sim.rda")
