library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("../sim/full_param2.R")

load("../data/full_boot.rda")

true.gen <- 1/sigma + 1/gamma
true.R <- beta/gamma
true.r <- 1/2 *(-(sigma+gamma)+sqrt((sigma-gamma)^2+4*beta*sigma))

true.par <- data.frame(
    type=c("gen", "r", "R"),
    value=c(true.gen, true.r, true.R)
)

rdf <- rlist %>%
    lapply(bind_rows, .id="run") %>%
    bind_rows(.id="sim")

coverage <- rdf %>%
    group_by(type) %>%
    summarize(coverage=mean(coverage))

rcover <- rdf %>%
    filter(type=="r", coverage==TRUE) %>%
    select(sim) %>%
    unlist %>%
    as.numeric

coverage2 <- rdf %>%
    filter(sim %in% rcover) %>%
    group_by(type) %>%
    summarize(coverage=mean(coverage))

g1 <- ggplot(coverage, aes(type, coverage)) +
    geom_point() +
    geom_hline(yintercept=0.95, lty=2) +
    scale_y_continuous(limits=c(0, 1)) +
    theme(axis.title.x = element_blank())

g2 <- g1 %+% (coverage2 %>% filter(type != "r"))

ggplot(rdf) +
    geom_boxplot(aes(type, estimate), width=0.3) +
    geom_hline(data = true.par, aes(yintercept=value), lty=2) +
    facet_wrap(~type, scale="free") +
    theme(axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

ggplot(rdf, aes(estimate, estimate/estimate)) +
    geom_point() +
    geom_errorbar(aes(ymin=lwr/estimate, ymax=upr/estimate)) +
    facet_wrap(~type ,scale="free")

ggplot(rdf, aes(estimate, coverage)) + 
    geom_point() +
    facet_wrap(~type ,scale="free")
