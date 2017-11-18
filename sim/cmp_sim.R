library(igraph)
source("cmp_param.R")
source("../R/seir.R")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

nsim <- 50

set.seed(101)
reslist <- vector('list', nsim)
i <- 1
while (i <= nsim) {
    print(i)
    rr <- seir(cmpGraph, beta, sigma, gamma, I0 = 10)
    if(nrow(rr$data) > 1000) {
        reslist[[i]] <- rr
        i <- i+1
    }
}

save("reslist", file="cmp_sim.rda")
