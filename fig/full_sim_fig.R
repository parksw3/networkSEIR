library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../sim/full_param.R")
source("../R/generation.R")
source("../R/empirical.R")

load("../data/full_sim.rda")

censor.gi <- lapply(
    reslist
    , network.generation
    , plot=FALSE
    , interval=c(1, 600)
    , interval.type="cases"
)

r <- lapply(
    reslist
    , function(x) {
        data <- x$data
        fdata <- data %>% filter(infected <= 600, infected >= 200)
        ll <- lm(log(infected)~time, data=fdata) 
        data.frame(r=ll$coefficients[2])
    }
)

true.r <- 1/2 *(-(sigma+gamma)+sqrt((sigma-gamma)^2+4*beta*sigma))

generation <- (
    censor.gi 
    %>% lapply(function(x) data.frame(interval=x)) 
    %>% bind_rows(.id="sim") 
    %>% merge(r %>% bind_rows(.id="sim")) 
    %>% as.tbl 
    %>% group_by(sim)
    %>% mutate(weight=exp(r*interval))
    %>% mutate(weight=weight/sum(weight))
)

intrinsic_fun <- function(tau) {
    sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

observed_fun <- function(tau) {
    beta*sigma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau)) * exp(-true.r*tau)
}

empty.df <- data.frame(
    x=c(0.1, 0.1)
    , y=c(0.1, 0.1)
    , group=c("intrinsic", "observed")
)

gg_base <- (
    ggplot(generation)
    + scale_x_continuous(name="Generation interval", limits=c(0, 15)) 
    + scale_y_continuous(limits=c(0, 0.7))
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
        , col='black', fill='#e7298a'
        , alpha=0.5, boundary=0, bins=30) 
    + geom_line(data=empty.df, aes(x, y, col=group, lty=group), lwd=1.2)
    + stat_function(fun=observed_fun, lwd=1.2, lty=2, xlim=c(0,15), col="#e7298a")
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=1, xlim=c(0,15), col="#d95f02")
    + scale_color_manual(values=c( "#d95f02", "#e7298a"))
    + ggtitle("Observed GI distributions")
    + theme(
        legend.position = c(0.85, 0.85)
        , legend.title = element_blank()
    )
)

gg2 <- (
    gg_base
    + geom_histogram(
        aes(interval, y=..density.., weight=weight)
        , col='black'
        , alpha=0.5, boundary=0, bins=30, fill="#d95f02")
    + stat_function(fun=observed_fun, lwd=1.2, lty=2, xlim=c(0,15), col="#e7298a")
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=1, xlim=c(0,15), col="#d95f02")
    + ggtitle("Corrected GI distributions")
)

gg_correction <- arrangeGrob(gg1, gg2, nrow=1)

ggsave("full_corrected_GI.pdf", gg_correction, width=6, height=3)
