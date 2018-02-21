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
    bind_rows(.id="sim") %>%
    mutate(key=factor(key, labels=c("corrected", "uncorrected"))) %>%
    filter(!(type=="r" & key=="uncorrected"))

coverage <- rdf %>%
    group_by(type, key) %>%
    summarize(coverage=mean(coverage)) %>%
    mutate(
        lwr=coverage-1.96*sqrt(coverage * (1-coverage)/100),
        upr=coverage+1.96*sqrt(coverage * (1-coverage)/100)
    )

rcover <- rdf %>%
    filter(type=="r", coverage==TRUE) %>%
    select(sim) %>%
    unlist %>%
    as.numeric

coverage2 <- rdf %>%
    filter(sim %in% rcover, key=="corrected") %>%
    group_by(type, key) %>%
    summarize(coverage=mean(coverage)) %>%
    mutate(
        lwr=coverage-1.96*sqrt(coverage * (1-coverage)/100),
        upr=coverage+1.96*sqrt(coverage * (1-coverage)/100)
    )

g1.base <- ggplot(NULL, aes(type, coverage)) +
    geom_hline(yintercept=0.95, lty=2) +
    scale_y_continuous(limits=c(0, 1)) +
    theme(axis.title.x = element_blank(),
          legend.title = element_blank(),
          legend.position = c(0.3, 0.8))

g1 <- g1.base %+% (coverage %>% filter(type != "r")) +
    geom_point(aes(col=key), size=2, position=position_dodge(width=0.2)) +
    geom_errorbar(aes(ymin=lwr, ymax=upr, col=key), width=0.2, position=position_dodge(width=0.2))
    

g1.r <- g1.base %+% (coverage %>% filter(type == "r")) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=lwr, ymax=upr), width=0.2)

gcover <- arrangeGrob(g1, g1.r, widths=c(0.4, 0.2))

ggsave("coverage.pdf", gcover, width=6, height=4)
  
gcover2 <- g1 %+% (coverage2 %>% filter(type != "r")) +
    theme(legend.position="none")

ggsave("coverage2.pdf", gcover2, width=6, height=4)

gest <- ggplot(rdf) +
    geom_boxplot(aes(type, estimate, col=key), width=0.3) +
    geom_hline(data = true.par, aes(yintercept=value), lty=2) +
    facet_wrap(~type, scale="free") +
    theme(axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())

ggsave("estimate.pdf", gest, width=6, height=4)

