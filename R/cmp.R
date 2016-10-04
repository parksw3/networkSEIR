library(igraph)
source("model.R")
cmpGraph <- read.table("../data/CA-CondMat.txt") ##importing condensed matter physics data
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)

pars <- c(
    beta = 0.01,
    sigma = 1/2,
    gamma = 1
)

r <- seir.gillespie(cmpGraph, pars, verbose = TRUE)
