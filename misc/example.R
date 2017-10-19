library(igraph)
source("../R/heap.R")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

## g <-sample_pa(n=500, power=1, m=20,  directed=F)

## randomly selected beta and gamma values. 
## Note that beta here is the contact rate between two individuals. 
beta <- 0.3
gamma <- 3

v_dist <- degree(cmpGraph)

mu_v <- mean(v_dist)

kappa <- var(v_dist)/mu_v + mu_v - 1

system.time(res <- seir.heap(cmpGraph, beta, gamma, seed=101))

## what the epidemic looks like
plot(res$data, type="l")

## calculated by considering *all* contacts
plot(density(unlist(res$intrinsic_generation)), main="intrinsic generation") 
curve(gamma*exp(-gamma*x), col=3, add=TRUE) ## looks pretty good

## calculated by considering only the first contacts
plot(density(res$t_infected-res$t_infected[res$infected_by], na.rm=TRUE), main="spatial generation") 
lines(density(unlist(res$intrinsic_generation)), main="intrinsic generation", col=2, lty=2) 
legend(1.5, 2.5, c("spatial", "intrinsic"), col=c(1, 2), lty=c(1, 2))
## It's not so clear what the spatial generation equation should be... I need to think about that a bit more
## I feel like Trapman doesn't quite work and we need to account for the fact that individuals with more neighbours
## will have higher rate of contact in this model
## we can at least see that the spatial generation is shorter than the intrinsic

## just checking that the infectious period from the simulation matches the expected distribution
plot(density(res$t_recovered-res$t_infected, na.rm=TRUE), col=2, main="infectious period")
curve(gamma*exp(-gamma*x), col=3, add=TRUE)
