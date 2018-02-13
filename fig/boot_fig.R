library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())

## why can't i push

load("../data/full_boot.rda")

rdf <- rlist %>%
    lapply(bind_rows, .id="run") %>%
    bind_rows(.id="sim")

coverage <- rdf %>%
    group_by(type) %>%
    summarize(coverage=mean(coverage))

ggplot(coverage, aes(type, coverage)) +
    geom_point() +
    geom_hline(yintercept=0.95, lty=2) +
    scale_y_continuous(limits=c(0, 1)) +
    theme(axis.title.x = element_blank())

ggplot(rdf, aes(type, estimate)) +
    geom_boxplot(width=0.3) +
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

    
    