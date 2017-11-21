library(igraph)

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

set.seed(101)
g <- sample_degseq(degree(cmpGraph), method = c("simple.no.multiple"))

save("g", file="random_network.rda")
