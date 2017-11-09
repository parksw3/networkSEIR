library(igraph)
source("../R/heap.R")
source("../R/generation.R")

branch <- 3
g <- graph.tree(branch^10, branch)

g <- graph.full(1000)

## randomly selected beta and gamma values. 
## Note that beta here is the contact rate between two individuals. 
beta <- 0.002
sigma <- Inf
gamma <- 1

system.time(res <- seir(g, beta, sigma, gamma, initial_infected = 1, seed=105))

## what the epidemic looks like
plot(res$data, type="l")

## calculated by considering *all* contacts
intrinsic.generation(res, breaks=20)
curve(gamma*exp(-gamma*x), add=TRUE) ## looks pretty good

di <- 0.2
ii <- seq(0, 6.5, by=di/2)

forward.mean <- backward.mean <- censor.forward.mean <- censor.backward.mean <- rep(NA, length(ii))

for (i in 1:length(ii)) {
    cf <- network.generation(res, plot=FALSE, type="forward", interval=c(0, ii[i]+di))
    cb <- network.generation(res, plot=FALSE, type="backward", interval=c(0, ii[i]+di))
    f <- network.generation(res, plot=FALSE, type="forward", interval=c(ii[i], ii[i]+di))
    b <- network.generation(res, plot=FALSE, type="backward", interval=c(ii[i], ii[i]+di))
    censor.forward.mean[i] <- mean(cf)
    censor.backward.mean[i] <- mean(cb)
    forward.mean[i] <- mean(f)
    backward.mean[i] <- mean(b)
}

pdf("interval_comparison.pdf", height=6, width=8)
par(mfrow=c(2, 2))
par(mar=c(4.2, 4.2, 4.2, 4.2))

plot(ii, censor.forward.mean, main="Forward", xlab="time", ylab="Censor based")
abline(h=1/gamma)
abline(h=1/(gamma+beta), lty=2)

plot(ii, censor.backward.mean, main="Backward", xlab="time", ylab=NA)
abline(h=1/gamma)
abline(h=1/(gamma+beta), lty=2)

plot(ii, forward.mean, xlab="time", ylab="Cohort based")
abline(h=1/gamma)
abline(h=1/(gamma+beta), lty=2)

plot(ii, backward.mean, xlab="time", ylab=NA)
abline(h=1/gamma)
abline(h=1/(gamma+beta), lty=2)

dev.off()

network.generation(res, breaks=20)
curve(gamma*exp(-gamma*x), add=TRUE)
curve((beta+gamma)*exp(-(beta+gamma)*x), add=TRUE, lty=2, col=2)
