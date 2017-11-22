library(igraph)
source("cmp_param.R")
source("../R/seir.R")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

d <- degree(cmpGraph)

reslist <- vector('list', length(d))
names(reslist) <- d

set.seed(101)
for (i in 1:length(d)) {
    print(i)
    rr <- seir.local(d[i]-1, beta, sigma, gamma)
    reslist[[i]] <- rr
}

save("reslist", file="local_sim.rda")
