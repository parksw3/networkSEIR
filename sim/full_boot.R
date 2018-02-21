library(bbmle)
library(MASS)
library(tidyr)
library(dplyr)

source("../sim/full_param2.R")
source("../R/generation.R")
source("../R/empirical.R")
source("../R/boot.R")

load("../data/full_sim2.rda")

verbose <- FALSE

true.gen <- 1/sigma + 1/gamma
true.R <- beta/gamma
true.r <- 1/2 *(-(sigma+gamma)+sqrt((sigma-gamma)^2+4*beta*sigma))

rlist <- vector('list', 100)

set.seed(101)
for (i in 1:100) {
    cat(i)
    sim <- reslist[[i]]
    rlist[[i]] <- replicate(1, {
        data <- generation.data(sim, imin=20, imax=200)
        bres <- generation.bootstrap(data)
        
        bres$coverage <- c(
            bres[1,3] < true.gen && true.gen < bres[1,4]
            , bres[2,3] < true.r && true.r < bres[2,4]
            , bres[3,3] < true.R && true.R < bres[3,4]
            , bres[4,3] < true.gen && true.gen < bres[4,4]
            , bres[5,3] < true.r && true.r < bres[5,4]
            , bres[6,3] < true.R && true.R < bres[6,4]
        )
        
        if (verbose) print(bres)
            
        bres
    }, simplify=FALSE)
}

save("rlist", file="full_boot.rda")
