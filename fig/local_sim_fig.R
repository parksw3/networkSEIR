library(bbmle)
library(igraph)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/cmp_param.R")

load("../data/local_sim.rda")

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

d <- degree(cmpGraph)

kappa <- var(d)/mean(d) + mean(d) -1

generation <- data.frame(
    interval=unlist(reslist)
)

intrinsic_fun <- function(tau) {
    sigma*(gamma)/(gamma-sigma) * (exp(-sigma*tau)-exp(-(gamma)*tau))
}

network_fun <- function(tau) {
    sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau))
}

R0_fun <- function(r) {
    (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)
}

empty.df <- data.frame(
    x=c(0.1, 0.1)
    , y=c(0.1, 0.1)
    , group=factor(c("intrinsic", "local"), level=c("intrinsic", "local"))
)

gg_local <- (
    ggplot(generation)
    + geom_histogram(
        aes(interval, y=..density..)
        , col='black', fill="#d95f02"
        , alpha=0.5, boundary=0, bins=30) 
    + geom_line(data=empty.df, aes(x, y, col=group), lwd=1.2)
    + ggtitle("Network GI distributions")
    + stat_function(fun=network_fun, lwd=1.2, lty=1, xlim=c(0,15), col="#d95f02")
    + stat_function(fun=intrinsic_fun, lwd=1.2, xlim=c(0,15), col="#66a61e")
    + scale_x_continuous(name="Generation interval", limits=c(0, 15)) 
    + scale_y_continuous(limits=c(0, 0.7))
    + scale_color_manual(values=c("#66a61e", "#d95f02"))
    + theme(
        panel.grid = element_blank()
        , panel.border=element_blank()
        , axis.line.x=element_line()
        , axis.title.y = element_blank()
        , axis.text.y = element_blank()
        , axis.ticks.y = element_blank()
        , plot.title = element_text(hjust=0.5)
        , legend.position = c(0.75, 0.85)
        , legend.title = element_blank()
    )
)

r <- seq(0, 3, by=0.05)

link <- data.frame(
    r=r,
    intrinsic=sapply(r, function(x) 1/integrate(function(y) intrinsic_fun(y) * exp(-x*y), 0, Inf)[[1]]),
    `spatially corrected`=R0_fun(r)
) %>%
    gather(key, value, -r)

gg_R <- ggplot(link) +
    geom_line(aes(r, value, col=key), lwd=1.2) +
    scale_x_continuous(name="Exponential growth rate", expand=c(0,0)) +
    scale_y_continuous(name="Reproductive number", expand=c(0,0), limits=c(0,20)) +
    geom_segment(data=data.frame(x=2, 
                                 y=14.5,
                                 x1=2, 
                                 y1=10),
                 aes(x, y, xend=x1, yend=y1),
                 arrow = arrow(length = unit(0.03, "npc"), type="closed"), lwd=0.9) +
    ggtitle("R-r relationship") +
    scale_color_manual(values=c("#66a61e", "#d95f02")) +
    theme(
        panel.grid = element_blank()
        , panel.border=element_blank()
        , axis.line=element_line()
        , plot.title = element_text(hjust=0.5)
        , legend.position = "none"
    )

gg <- arrangeGrob(gg_local, gg_R, nrow=1)

ggsave("local.pdf", gg, width=6, height=3)
