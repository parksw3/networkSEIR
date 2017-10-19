library(igraph)
source("../R/heap.R")

## a scale-free graph
g <-sample_pa(n=500, power=1, m=20,  directed=F)

## randomly selected beta and gamma values. 
## Note that beta here is the contact rate between two individuals. 
## So the homogeneous (everyone is connected to everyon) equivalent beta would be beta * (number of verticies-1)  
## We need this homogeneous equivalent beta when calculating the spatial generation distribution
## Later, we want to fix beta for each individual and distribute contact rates among neighbours
beta <- 0.4
gamma <- 3

v_dist <- degree(g)

mu_v <- mean(v_dist)

kappa <- var(v_dist)/mu_v + mu_v - 1

res <- seir.heap(g, beta, gamma, seed=101)

## what the epidemic looks like
plot(res$data, type="l")

## calculated by considering *all* contacts
plot(density(unlist(res$forward_generation)), main="intrinsic generation") 
curve(gamma*exp(-gamma*x), col=3, add=TRUE) ## looks pretty good

## calculated by considering only the first contacts
plot(density(res$t_infected-res$t_infected[res$infected_by], na.rm=TRUE), main="spatial generation") 
curve(gamma*exp(-gamma*x), col=2, add=TRUE)
curve((beta*(length(V(g))-1)/kappa+gamma)*exp(-(beta*(length(V(g))-1)/kappa+gamma)*x), col=3, add=TRUE)
## I think homogeneous equivalent beta has to appear here
## I would say that this looks pretty good

## just checking that the infectious period from the simulation matches the expected distribution
plot(density(res$t_recovered-res$t_infected, na.rm=TRUE), col=2, main="infectious period")
curve(gamma*exp(-gamma*x), col=3, add=TRUE)
