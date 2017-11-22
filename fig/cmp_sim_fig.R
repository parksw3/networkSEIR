library(bbmle)
library(igraph)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/cmp_param.R")
source("../R/generation.R")
source("../R/empirical.R")
source("../R/mle.R")

load("../data/cmp_sim.rda")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

d <- degree(cmpGraph)

kappa <- var(d)/mean(d) + mean(d) -1

censor.gi <- lapply(
    reslist
    , network.generation
    , plot=FALSE
    , interval=c(200, 400)
    , interval.type="cases"
)

empirical <- (
    lapply(reslist, empirical.R0, n=200)
    %>% lapply(function(x) data.frame(empirical=x))
    %>% bind_rows(.id="sim")
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

generation <- (
    censor.gi 
    %>% lapply(function(x) data.frame(interval=x)) 
    %>% bind_rows(.id="sim") 
    %>% merge(r %>% bind_rows(.id="sim")) 
    %>% as.tbl 
    %>% group_by(sim)
    %>% mutate(weight=exp(r*interval)/sum(exp(r*interval)))
)

R0 <- (
    generation 
    %>% group_by(sim) 
    %>% summarize(
        r=mean(r)
        , uncorrected=1/mean(exp(-r*interval))
        , corrected=mean(exp(r*interval)))
    %>% merge(empirical)
    %>% mutate(intrinsic=R0fun(r),
               network=R0fun.network(r))
    %>% gather(key, value, -sim, -r)
    %>% mutate(key=factor(key, levels=c("uncorrected", "corrected", "empirical", "network", "intrinsic")))
    %>% as.tbl
)

intrinsic_fun <- function(tau) {
    sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau))
}

observed_fun <- function(tau) {
    true.R0*sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau)) * exp(-mean.r*tau)
}

empty.df <- data.frame(
    x=c(0.1, 0.1)
    , y=c(0.1, 0.1)
    , group=c("censored (network)", "intrinsic (network)")
)
    
gg_base <- (
    ggplot(generation)
    + scale_x_continuous(name="time") 
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

gg1 <- (
    gg_base
    + geom_histogram(
        aes(interval, y=..density..)
        , col='black', fill='grey'
        , alpha=0.7, boundary=0, bins=30) 
    + geom_line(data=empty.df, aes(x, y, lty=group), lwd=1.2)
    + stat_function(fun=observed_fun, lwd=1.2, xlim=c(0,10))
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=2, alpha=0.2, xlim=c(0,10))
    + ggtitle("Observed GI distributions during exponential growth phase")
    + theme(
        legend.position = c(0.85, 0.85)
        , legend.title = element_blank()
    )
)

gg2 <- (
    gg_base
    + geom_histogram(
        aes(interval, y=..density..)
        , col='grey', fill='black'
        , alpha=0.15, boundary=0, bins=30) 
    + stat_function(fun=observed_fun, lwd=1.2, alpha=0.1, xlim=c(0,10))
    + geom_histogram(
        aes(interval, y=..density.., weight=weight)
        , fill='grey', col='black'
        , alpha=0.7, boundary=0, bins=30)
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=2, xlim=c(0,10))
    + ggtitle("Corrected GI distributions")
)

gg_R <- (
    ggplot(R0, aes(key, value)) 
    + geom_boxplot(alpha=0.5, width=0.7)
    + scale_y_log10("Reproductive number", breaks=c(3, 5, 10, 20))
    + theme(
        panel.grid=element_blank(),
        axis.title.x=element_blank()
    )
)

gg_correction <- arrangeGrob(gg1, gg2, gg_R, layout_matrix=cbind(c(1,2), c(3,3)), widths=c(0.6, 0.4))

ggsave("corrected_GI.pdf", gg_correction, width=10, height=6)
