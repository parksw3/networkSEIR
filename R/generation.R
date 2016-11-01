## Recreate figure 5 from Trapman et al paper
source("parameters.R")
source("functions.R")
load("../data/condmat_sim.rda")

l <- lapply(datList, function(x){
    gen <- generation(x, interval = c(270, 300))
})

dist <- unlist(l)

plot(density(dist))
with(as.list(pars),{
    curve(sigma * (beta + gamma)/(beta + gamma - sigma) *
              (exp(- sigma * x) - exp(-(beta+gamma)*x)), add = TRUE,
          lty = 2)
})

with(as.list(pars),{
    curve(sigma * gamma/(gamma - sigma) *
              (exp(- sigma * x) - exp(-(gamma)*x)), add = TRUE,
          lty = 3)
})
