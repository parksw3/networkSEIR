library(igraph)
library(data.table)
source("model.R")
cmpGraph <- read.table("../data/CA-CondMat.txt") ##importing condensed matter physics data
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

## this parameter gives R0 of approximately 2
## These values seem to match orange vertical line of figure 5 

pars <- c(
    beta = 0.02,
    sigma = 1/9.4,
    gamma = 1/5,
    I0 = 10
)

## example simulation using this parameter

fn <- "cmpSim1.rda"

if(!file.exists(fn)){
    cmpGraph.sim <- seir.gillespie(cmpGraph, pars, verbose = TRUE,
                                   seed = 101)
    
    save("cmpGraph.sim", file = fn)
}else{
    load(fn)
}

## 10 simulations due to lack of computation horsepower...

fn2 <- "cmpSim10.rda"

if(!file.exists(fn2)){
    set.seed(101)
    sumList <- list()
    datList <- list()
    i <- 1
    
    while(i <= 10){
        cat(i)
        sim <- seir.gillespie(cmpGraph, pars)
        if(!is.na(sim)){
            sumList[[i]] <- sim$epidemic.summary
            print(sumList[[i]])
            datList[[i]] <- sim$epidemic.data
            i <- i + 1
        }
    }
    save("sumList", "datList", file = fn2)
}else{
    load(fn2)
}
