library(igraph)
source("model.R")
cmpGraph <- read.table("../data/CA-CondMat.txt") ##importing condensed matter physics data
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

pars <- c(
    beta = 0.15,
    sigma = 1/2,
    gamma = 1,
    I0 = 10
)

set.seed(101)
cmpGraph.sim <- seir.gillespie(cmpGraph, pars, verbose = TRUE)

epi <- cmpGraph.sim$result
matplot(epi[,1], epi[,-1], type = "l")

gen <- cmpGraph.sim$generation
hist(gen, breaks = c(0:25))
mean(gen)
var(gen)

set.seed(101)
r <- seir.gillespie(g, pars)
mean(r$generation)
var(r$generation)
hist(r$generation, breaks = c(0:30))

r2 <- r$result

