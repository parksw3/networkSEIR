source("../sim/full_param.R")
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

source("../R/generation.R")
load("../data/full_sim.rda")

censor.gi <- lapply(
    reslist
    , network.generation
    , plot=FALSE
    , interval=c(2, 500)
    , interval.type="cases"
)

r <- lapply(
    reslist
    , function(x) {
        data <- x$data
        fdata <- data %>% filter(infected <= 500, infected >= 250)
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
    %>% mutate(weight=exp(r*interval))
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
    , group=c("censored", "intrinsic")
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
    + stat_function(fun=observed_fun, lwd=1.2, xlim=c(0,11))
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=2, alpha=0.2, xlim=c(0,11))
    + ggtitle("Observed GI distributions during exponential growth phase")
    + theme(
        legend.position = c(0.85, 0.85)
        , legend.title = element_blank()
    )
)

gg2 <- (
    gg_base
    + stat_function(fun=observed_fun, lwd=1.2, alpha=0.1, xlim=c(0,11))
    + geom_histogram(
        aes(interval, y=..density.., weight=weight)
        , fill='grey', col='black'
        , alpha=0.7, boundary=0, bins=40)
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=2, xlim=c(0,11))
    + ggtitle("Corrected GI distributions")
)

R0 <- (
    generation 
    %>% group_by(sim) 
    %>% summarize(uncorrected=1/mean(exp(-r*interval)),
        corrected=mean(exp(r*interval))) 
    %>% gather(key, value, -sim)
    %>% mutate(key=factor(key, level=c("uncorrected", "corrected")))
)
    
gg_R <- (
    ggplot(R0) 
    + geom_boxplot(aes(key, value), fill="grey", alpha=0.5)
    + scale_y_continuous("Reproductive number")
    + geom_hline(yintercept=beta/gamma, lty=2)
    + theme(
        panel.grid=element_blank(),
        axis.title.x=element_blank()
    )
)

gg_correction <- arrangeGrob(gg1, gg2, gg_R, layout_matrix=cbind(c(1,2), c(3,3)), widths=c(1, 0.4))

ggsave("corrected_GI.pdf", gg_correction, width=8, height=6)
