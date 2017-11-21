library(igraph)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(grid)
library(gridExtra)

source("../sim/cmp_param.R")
source("../R/generation.R")
source("../R/empirical.R")
source("../R/mle.R")

load("../data/random_sim.rda")
load("../data/random_network.rda")

d <- degree(g)

kappa <- var(d)/mean(d) + mean(d) -1

censor.gi <- lapply(
    reslist
    , network.generation
    , plot=FALSE
    , interval=c(200, 400)
    , interval.type="cases"
)

r <- lapply(
    reslist
    , function(x) {
        data <- x$data
        fdata <- data %>% filter(infected <= 400, infected >= 200)
        ll <- lm(log(infected)~time, data=fdata) 
        data.frame(r=ll$coefficients[2])
    }
)

mean.r <- mean(unlist(r))

R0fun <- function(r) (1+r/gamma)*(1+r/sigma)

R0fun.network <- function(r) (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)

intrinsic.R0 <- R0fun(mean.r)

true.R0 <- R0fun.network(mean.r)

quant.r <- quantile(unlist(r), c(0.05, 0.95))

generation <- (
    censor.gi 
    %>% lapply(function(x) data.frame(interval=x)) 
    %>% bind_rows(.id="sim") 
    %>% merge(r %>% bind_rows(.id="sim")) 
    %>% as.tbl 
    %>% group_by(sim)
    %>% mutate(weight=exp(r*interval))
    %>% mutate(weight=weight/sum(weight))
    %>% filter(r > quant.r[1], r < quant.r[2])
)

intrinsic_fun <- function(tau) {
    sigma*(gamma)/(gamma-sigma) * (exp(-sigma*tau)-exp(-(gamma)*tau))
}

network_fun <- function(tau) {
    sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau))
}

network_observed_fun <- function(tau) {
    true.R0*sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau)) * exp(-mean.r*tau)
}

gg_base <- (
    ggplot(generation)
    + scale_x_continuous(name="time") 
    + scale_y_continuous(limits=c(0, 1))
    + theme(
        panel.grid = element_blank()
        , panel.border=element_blank()
        , axis.line.x=element_line()
        , axis.title.y = element_blank()
        , axis.text.y = element_blank()
        , axis.ticks.y = element_blank()
        , plot.title = element_text(hjust=0.5)
    )
)

g_intrinsic <- gg_base + 
    geom_histogram(aes(interval, y=..density.., weight=weight), 
                   col='black', fill='grey', 
                   alpha=0.7, boundary=0, bins=15) +
    geom_vline(xintercept=1/sigma+1/gamma, lwd=1) +
    geom_vline(xintercept=weighted.mean(generation$interval, generation$weight), lwd=1, lty=2) +
    stat_function(fun=intrinsic_fun, lwd=1.2, xlim=c(0,8)) +
    ggtitle("intrinsic GI distributions")
    
g_network <- gg_base +
    geom_histogram(aes(interval, y=..density.., weight=weight), 
                   col='black', fill='grey', 
                   alpha=0.7, boundary=0, bins=15) +
    geom_vline(xintercept=1/sigma+1/(gamma+beta), lwd=1) +
    geom_vline(xintercept=weighted.mean(generation$interval, generation$weight), lwd=1, lty=2) +
    stat_function(fun=network_fun, lwd=1.2, xlim=c(0,8)) +
    ggtitle("network GI distributions")
    
g_observed <- gg_base +
    geom_histogram(aes(interval, y=..density..), 
                   col='black', fill='grey', 
                   alpha=0.7, boundary=0, bins=15) +
    geom_vline(xintercept = integrate(function(x) x*network_observed_fun(x), 0, Inf)[[1]], lwd=1) +
    geom_vline(xintercept=mean(generation$interval), lwd=1, lty=2) +
    stat_function(fun=network_observed_fun, lwd=1.2, xlim=c(0,8)) +
    ggtitle("observed GI distributions")

gg <- arrangeGrob(g_observed, g_network, g_intrinsic, nrow=1)

ggsave("random_sim.pdf", gg, width=8, height=2)
