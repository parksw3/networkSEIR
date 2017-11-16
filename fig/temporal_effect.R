source("../sim/full_param.R")
library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../R/generation.R")
load("../data/full_sim.rda")

## simulation

censor.t <- seq(0, 70, by=2)

censor.gi <- lapply(censor.t, function(t){
    lapply(
        reslist
        , network.generation
        , plot=FALSE
        , interval=c(0, t)
    )
})

names(censor.gi) <- censor.t

r <- lapply(
    reslist
    , function(x) {
        data <- x$data
        fdata <- data %>% filter(time <= 20, time >= 10)
        ll <- lm(log(infected)~time, data=fdata) 
        data.frame(r=ll$coefficients[2])
    }
)

true.r <- 1/2 *(-(sigma+gamma)+sqrt((sigma-gamma)^2+4*beta*sigma))

generation <- censor.gi %>%
    lapply(lapply, function(x) data.frame(interval=x)) %>%
    lapply(bind_rows, .id="sim") %>% 
    bind_rows(.id="t") %>%
    as.tbl

mean.generation <- generation %>%
    group_by(t) %>%
    summarize(mean=mean(interval)) %>%
    mutate(t=as.numeric(t))

exp.generation <- generation %>%
    filter(t==20) %>%
    merge(r %>% bind_rows(.id="sim")) %>%
    as.tbl

intrinsic_fun <- function(tau) {
    sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

observed_fun <- function(tau) {
    beta*sigma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau)) * exp(-true.r*tau)
}

## ode

seir.ode <- function(t, y, parms) {
    with(as.list(c(y, parms)), {
        dS <- -beta/N*S*I
        dE <- beta/N*S*I - sigma*E
        dI <- sigma*E - gamma*I
        list(c(dS, dE, dI))
    })
}

parms <- c(beta=beta, sigma=sigma, gamma=gamma, N=N)
y <- c(S = N-1, E=0, I=1)

dt <- 0.05
tmax <- 70
t <- seq(0, tmax, dt)
out <- as.data.frame(ode(y, t, seir.ode, parms))
names(out)[1] <- "t"

incidence <- beta*out$S*out$I/N
c.incidence <- cumsum(incidence)

out$incidence <- incidence

gen <- intrinsic_fun(t)

censor <- function(tau, tmax) {
    m <- which(t==tmax)
    total <- c.incidence[m]
    
    i <- which(t==tau)
    
    beta*sum(incidence[1:(m-i+1)] * out$S[i:m]/N)/total*gen[i]
}

mean.censor <- function(tmax) {
    if (tmax<=0){
        return(0)
    } else{
        t <- t[t <= tmax]
        
        dc <- sapply(t, censor, tmax=tmax)
        return(sum(dc*t)/sum(dc))
    }
}

tau <- seq(0, tmax, 0.5)
cm <- sapply(tau, mean.censor)

censor.df <- data.frame(
    t=tau
    , mean=cm
)

exp.df <- data.frame(
    t=t[which(t==20)],
    incidence=incidence[which(t==20)],
    mean=cm[which(tau==20)]
)

out.shade <- rbind(out[out$t <= 20,c("t", "incidence")], c(20, 0), c(0,0))

g1 <- ggplot(out, aes(t, incidence)) +
    geom_line(lwd=1.2) +
    geom_polygon(data=out.shade, alpha=0.2) +
    geom_segment(data=data.frame(t=20, incidence=0, tend=20, iend=out$incidence[which(t==20)]), aes(xend=tend, yend=iend),lty=2) +
    geom_point(data=exp.df, size=4) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0), limits=c(0, 1800)) +
    theme(
        panel.grid=element_blank(),
        axis.title.x=element_blank()
    )

g2 <- ggplot(censor.df, aes(t, mean)) +
    geom_line(lwd=1.2) +
    geom_point(data=exp.df, size=4) +
    geom_point(data=mean.generation, shape=1, size=2) +
    geom_segment(data=data.frame(t=20, mean=0, tend=20, iend=censor.df$mean[which(censor.df$t==20)]), aes(xend=tend, yend=iend),lty=2) +
    geom_hline(yintercept = 1/sigma + 1/gamma, lty=2) +
    scale_x_continuous("time", expand=c(0,0)) +
    scale_y_continuous("mean censored GI", expand=c(0,0), limits=c(0, 4)) +
    theme(
        panel.grid=element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))
    )

empty.df <- data.frame(
    x=c(0.1, 0.1)
    , y=c(0.1, 0.1)
    , group=c("censored", "intrinsic")
)

g3 <- ggplot(exp.generation) +
    geom_histogram(aes(interval, y=..density..), col='black', fill='grey', alpha=0.8, boundary=0, bins=40) + 
    geom_line(data=empty.df, aes(x, y, lty=group), lwd=1.2) +
    xlab("generation interval") +
    stat_function(fun=observed_fun, lwd=1.2, xlim=c(0,14)) +
    stat_function(fun=intrinsic_fun, lwd=1.2, lty=2, xlim=c(0,14)) +
    ggtitle("Observed GI distributions at t=20") +
    theme(
        panel.grid = element_blank()
        , panel.border=element_blank()
        , axis.line.x=element_line()
        , axis.title.y = element_blank()
        , axis.text.y = element_blank()
        , axis.ticks.y = element_blank()
        , plot.title = element_text(hjust=0.5)
        , legend.position = c(0.85, 0.9)
        , legend.title = element_blank()
    )

gg_temporal <- arrangeGrob(g1, g2, g3, layout_matrix=cbind(c(1,2), c(3,3)), widths=c(0.8, 0.6))

ggsave("temporal_effect.pdf", gg_temporal, width=8, height=6)

