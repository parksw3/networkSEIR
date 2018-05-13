source("ebola_param.R")
true.R <- beta/gamma

true.mean <- 1/sigma + 1/gamma

set.seed(101)
nQuant <- 10000
q <- (2*(1:nQuant)-1)/(2*nQuant)

lat <- qexp(q, rate=sigma)
inf <- qexp(q, rate=gamma)

ii <- sample(inf, prob=inf, replace=TRUE)

gen <- lat + runif(nQuant, min=0, max=ii)

## approximately true...
true.shape <- true.mean^2/var(gen)

gi.fun <- function(n) rgamma(n, shape=true.shape, scale=true.mean/true.shape)

nsim <- 100

set.seed(101)
reslist <- vector('list', nsim)
j <- 1
while (j <= nsim) {
    print(j)
    
    t_infected <- rep(0, 10)
    t_gen <- rep(NA, 10)
    i <- 1
    infected_by <- rep(NA, 10)
    
    while(i < 1000 && i <= length(t_infected)) {
        inf <- rpois(1, true.R)
        if (inf > 0) {
            gen <- gi.fun(inf)
            t_infected <- c(t_infected, gen + t_infected[i])
            infected_by <- c(infected_by, rep(i, inf))
            t_gen <- c(t_gen, gen)
            oo <- order(t_infected)
            t_infected <- t_infected[oo]
            t_gen <- t_gen[oo]
            infected_by <- infected_by[oo]
            
        }
        
        i <- i + 1
    }
    
    if (i==1000) {
        reslist[[j]] <- data.frame(
            index=1:1000,
            t_infected=t_infected[1:1000],
            infected_by=infected_by[1:1000],
            generation=t_gen[1:1000]
        )
        
        j <- j + 1
    }
}

save("reslist", file="poisson_sim.rda")
