source("../sim/full_param.R")
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

seir.ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
        dS <- -beta/N*S*I
        dE <- beta/N*S*I - sigma*E
        dI <- sigma*E - gamma*I
        list(c(dS, dE, dI))
    })
}

parms <- list(beta=2, sigma=0.5, gamma=1, N=20000)
y <- c(S = parms$N-1, E=0, I=1)

dt <- 0.01
tmax <- 80
t <- seq(0, tmax, dt)
out <- as.data.frame(ode(y, t, seir.ode, parms))
incidence <- parms$beta*out$S*out$I/parms$N*dt

gen <- function(tau) parms$sigma*parms$gamma/(parms$sigma - parms$gamma) * (exp(-parms$gamma*tau)-exp(-parms$sigma*tau))

cg <- function(tau, tmax) {
    t <- t[t<tmax]
    total <- sum(incidence)
    
    i <- which(t>tau)
    
    2*sum(incidence[1:length(i)] * out$S[i]/1000)/total*gen(tau)
}

ff <- function(tmax) {
    t <- t[t <= tmax]
    
    dc <- sapply(t, cg, tmax=tmax)
    sum(dc*t)/sum(dc)
}

tau <- seq(0.5, tmax, 0.5)
cm <- sapply(tau, ff)

plot(tau, cm, type="l", ylim=c(0,4), ylab="mean censored GI", xlab="time")
abline(h=1/parms$gamma+1/parms$sigma, lty=2)

