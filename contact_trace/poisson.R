R0 <- 2
gi.fun <- function(n) rgamma(n, shape=2, scale=4)

t.infected <- 0
t.gen <- 0
i <- 1

while(length(t.infected) < 10000 || i > length(t.infected)) {
    inf <- rpois(1, R0)
    if (inf > 0) {
        gen <- gi.fun(inf)
        t.infected <- c(t.infected, gen + t.infected[i])
        t.gen <- c(t.gen, gen)
        oo <- order(t.infected)
        t.infected <- t.infected[oo]
        t.gen <- t.gen[oo]
    }
    
    i <- i + 1
}

day <- round(t.infected)
plot(table(day))

hist(t.gen[t.infected < 70], freq=FALSE)

curve(dgamma(x, shape=2, scale=4), add=TRUE)
