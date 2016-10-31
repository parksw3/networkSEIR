library(igraph)
library(data.table)
source("model.R")
## importing condensed matter physics data
cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

## This parameter gives R0 of approximately 2
## These values seem to match orange vertical line of figure 5 

pars <- c(
    beta = 0.02,
    sigma = 1/9.4,
    gamma = 1/5,
    I0 = 10
)

fn <- "condmat_sim.rda"

if(!file.exists(fn)){
    set.seed(101)
    sumList <- list()
    datList <- list()
    i <- 1
    
    while(i <= 1000){
        cat(i)
        sim <- seir.gillespie(cmpGraph, pars)
        if(!is.na(sim)){
            sumList[[i]] <- sim$epidemic.summary
            print(sumList[[i]])
            datList[[i]] <- sim$epidemic.data
            i <- i + 1
        }
        save("sumList", "datList", file = fn)
    }
}else{
    load(fn)
}
