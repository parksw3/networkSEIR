library(igraph)
source("../R/heap.R")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

## cmpGraph <- sample_pa(n=1000, power=2, m=20,  directed=F)
## cmpGraph <- graph.full(1000)

## randomly selected beta and gamma values. 
## Note that beta here is the contact rate between two individuals. 
beta <- 0.1
sigma <- 0.5
gamma <- 1

v_dist <- degree(cmpGraph)

mu_v <- mean(v_dist)

kappa <- var(v_dist)/mu_v + mu_v - 1

system.time(res <- seir.heap(cmpGraph, beta, sigma, gamma, seed=101))

## what the epidemic looks like
plot(res$data, type="l")

## calculated by considering *all* contacts
hist(unlist(res$intrinsic_generation), main="intrinsic generation", freq=FALSE, breaks=50) 
curve(sigma*gamma/(sigma-gamma)*(exp(-gamma*x) - exp(-sigma*x)), add=TRUE) ## looks pretty good

## calculated by considering only the first contacts
hist(res$t_infected-res$t_infected[res$infected_by], main="spatial generation", freq=FALSE, breaks=50) 
curve(sigma*gamma/(sigma-gamma)*(exp(-gamma*x) - exp(-sigma*x)), add=TRUE) ## intrinsic
curve(sigma*(gamma+beta)/(sigma-(gamma+beta))*(exp(-(gamma+beta)*x) - exp(-sigma*x)), lty=2, add=TRUE) ## network

## just checking that the infectious period from the simulation matches the expected distribution
hist(res$t_recovered-res$t_infectious, main="infectious period", freq=FALSE, breaks=50)
curve(gamma*exp(-gamma*x), col=2, add=TRUE)

## latent period
hist(res$t_infectious-res$t_infected, main="infectious period", freq=FALSE, breaks=50)
curve(sigma*exp(-sigma*x), col=2, add=TRUE)
