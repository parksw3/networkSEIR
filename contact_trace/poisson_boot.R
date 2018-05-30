source("ebola_param.R")
source("../R/generation.R")
source("../R/empirical.R")
source("minuslogl.R")
source("bootfun.R")

load("poisson_sim.rda")

rlist <- vector('list', 3)

nsample <- c(100, 500)

set.seed(101)
for (j in 1:length(nsample)) {
    
    subrlist <- vector('list', 100)
    
    for (i in 1:100) {
        cat(i)
        
        data <- reslist[[i]][1:nsample[j]+10,]
        
        res <- bootfun(data)
        
        res$cover <- c(
            res$lwr < c(true.mean, true.shape, true.R) & c(true.mean, true.shape, true.R) < res$upr
        )
        
        rownames(res) <- NULL
        
        print(res)
        
        subrlist[[i]] <- res
    }
    
    rlist[[j]] <- subrlist
}

save("rlist", file="poisson_boot.rda")
