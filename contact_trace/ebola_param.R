beta <- 2/5
sigma <- 1/11.4
gamma <- 1/5
N <- 60000

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

true.R <- beta/gamma