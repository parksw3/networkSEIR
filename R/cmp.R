library(igraph)
source("model.R")
cmpGraph <- read.table("../data/CA-CondMat.txt") ##importing condensed matter physics data
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

pars <- c(
    beta = 0.03,
    sigma = 1/9,4,
    gamma = 1/5,
    I0 = 10
)

## set.seed(101)
## cmpGraph.sim <- seir.gillespie(cmpGraph, pars, verbose = TRUE)

## save("cmpGraph.sim", file = "cmpSim1.rda")

load("compSim1.rda")

print(cmpGraph.sim$summary)

epi <- cmpGraph.sim$result
matplot(epi[,1], epi[,"I"], type = "l")

gen <- cmpGraph.sim$generation
hist(gen, breaks = c(0:120))
mean(gen)
var(gen)

# Homogeneous case
#
# set.seed(101)
# r <- seir.gillespie(g, pars)
# mean(r$generation)
# var(r$generation)
# hist(r$generation, breaks = c(0:30))
# 
# r2 <- r$result
