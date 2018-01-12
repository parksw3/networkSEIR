source("../sim/full_param.R")
source("../R/generation.R")
source("../R/empirical.R")
source("../R/boot.R")

load("../data/full_sim.rda")

reporting.rate <- 1

true.gen <- 1/sigma+1/gamma

rlist <- vector('list', 100)

set.seed(101)
for (i in 1:100) {
    cat(i)
    sim <- reslist[[i]]
    rlist[[i]] <- replicate(5, {
        df <- sim
        index <- which(!is.na(sim$t_infected))
        reported <- runif(length(index)) < reporting.rate
        df$infected_by[index[!reported]] <- NA
        gen <- network.generation(df, plot=FALSE)
        
        df2 <- data.frame(
            time=sort(df$t_infected[index[reported]]),
            infected=1:sum(reported)
        )
        
        r <- coef(lm(log(infected)~time, data=df2))[2]
        
        weight <- exp(r*gen)
        
        b1 <- boot(gen, weight)
        b2 <- boot_weight(gen, weight)
        b3 <- boot_half(gen, weight)
        
        data.frame(
            name=c("random", "weighted", "half"),
            width=c(diff(b1), diff(b2), diff(b3)),
            coverage=c(b1[1]<true.gen&&true.gen<b1[2], b2[1]<true.gen&&true.gen<b2[2], b3[1]<true.gen&&true.gen<b3[2]),
            estimate=weighted.mean(gen, weight)
        )
    }, simplify=FALSE)
}

save("rlist", file="full_boot.rda")

