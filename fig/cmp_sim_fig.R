library(igraph)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/cmp_param.R")
source("../R/generation.R")
source("../R/reproductive.R")

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
    , interval=c(10, 200)
    , interval.type="cases"
)

empirical <- (
    lapply(reslist, empirical.R0, n=50)
    %>% lapply(function(x) data.frame(empirical=x))
    %>% bind_rows(.id="sim")
)

r <- lapply(
    reslist
    , function(x) {
        data <- x$data
        fdata <- data %>% filter(infected <= 200, infected >= 50)
        ll <- lm(log(infected)~time, data=fdata) 
        data.frame(r=ll$coefficients[2])
    }
)

findr <- function(r) (gamma+r)*(sigma+r)/((kappa-1)*sigma-r)-beta

true.r <- uniroot(findr, c(0, 4))$root

R0fun <- function(r) (1+r/gamma)*(1+r/sigma)

R0fun.network <- function(r) (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)

naive.R0 <- R0fun(true.r)

true.R0 <- R0fun.network(true.r)

quant.r <- quantile(unlist(r), c(0.05, 0.95))

generation <- (
    censor.gi 
    %>% lapply(function(x) data.frame(interval=x)) 
    %>% bind_rows(.id="sim") 
    %>% merge(r %>% bind_rows(.id="sim")) 
    %>% as.tbl 
    %>% group_by(sim)
    %>% mutate(weight=exp(r*interval)/sum(exp(r*interval)))
    %>% filter(r > quant.r[1], r < quant.r[2])
)

R0 <- (
    generation 
    %>% group_by(sim) 
    %>% summarize(
        r=mean(r)
        , uncorrected=1/mean(exp(-r*interval))
        , corrected=mean(exp(r*interval))) 
    %>% merge(empirical)
    %>% mutate(naive=R0fun(r))
    %>% gather(key, value, -sim, -r)
    %>% mutate(key=factor(key, levels=c("uncorrected", "corrected", "empirical", "naive")))
    %>% as.tbl
)

intrinsic_fun <- function(tau) {
    sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau))
}

observed_fun <- function(tau) {
    true.R0*sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau)) * exp(-true.r*tau)
}

empty.df <- data.frame(
    x=c(0.1, 0.1)
    , y=c(0.1, 0.1)
    , group=c("censored (network)", "intrinsic (network)")
)
    
gg_base <- (
    ggplot(generation)
    + geom_histogram(
        aes(interval, y=..density..)
        , col='grey', fill='black'
        , alpha=0.15, boundary=0, bins=40) 
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
    + stat_function(fun=observed_fun, lwd=1.2, alpha=0.1, xlim=c(0,10))
    + geom_histogram(
        aes(interval, y=..density.., weight=weight)
        , fill='grey', col='black'
        , alpha=0.7, boundary=0, bins=40)
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=2, xlim=c(0,10))
    + ggtitle("Corrected GI distributions")
)

gg_R <- (
    ggplot(R0, aes(key, value)) 
    + geom_boxplot(alpha=0.5, width=0.4)
    + scale_y_log10("Reproductive number", breaks=c(2, 5, 10, 50))
    + geom_hline(yintercept=naive.R0, lty=3)
    + geom_hline(yintercept=true.R0, lty=2)
    + theme(
        panel.grid=element_blank(),
        axis.title.x=element_blank()
    )
)

gg_correction <- arrangeGrob(gg1, gg2, gg_R, layout_matrix=cbind(c(1,2), c(3,3)), widths=c(0.8, 0.4))

ggsave("corrected_GI.pdf", gg_correction, width=10, height=6)
