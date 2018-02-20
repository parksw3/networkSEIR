source("full_param2.R")
source("../R/seir.R")
source("../R/generation.R")

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    rr <- seir.full(N, beta/N, sigma, gamma, I0 = 1, imax=10000)
    if(nrow(rr$data) > 8000) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="full_sim2.rda")
