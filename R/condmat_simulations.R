library(igraph)
library(data.table)
source("model.R")
source("parameters.R")
## importing condensed matter physics data
cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

fn <- "condmat_sim.rda"

set.seed(101)
    
n <- 200
    
sumList <- vector("list", n)
datList <- vector("list", n)
i <- 1

while(i <= n){
    cat(i)
    sim <- try(seir.gillespie(cmpGraph, pars))
    
    if(!inherits(sim, "try-error")){
        sumList[[i]] <- sim$epidemic.summary
        print(sumList[[i]])
        datList[[i]] <- sim$epidemic.data
        i <- i + 1
    }
    save("sumList", "datList", file = fn)
}