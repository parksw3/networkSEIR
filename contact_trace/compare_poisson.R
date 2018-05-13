library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("ebola_param.R")

mlelist <- vector('list', 2)

load("poisson_mle.rda")

mlelist[[1]] <- rlist

load("seir_mle.rda")

mlelist[[2]] <- rlist

df <- mlelist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="type") %>%
    mutate(type=factor(type, labels=c("gamma", "SEIR"))) %>%
    rename(likelihood=method)

g1 <- ggplot(df %>% filter(param=="R")) +
    geom_boxplot(aes(type, mean)) +
    geom_hline(yintercept=2, lty=2)+
    ylab("reproductive number") +
    xlab("")

g2 <- ggplot(df %>% filter(param=="mean")) +
    geom_boxplot(aes(type, mean, fill=likelihood), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.mean, lty=2) +
    ylab("mean generation (days)") +
    xlab("") +
    theme(
        legend.position = "none"
    )

g3 <- ggplot(df %>% filter(param=="shape")) +
    geom_boxplot(aes(type, 1/sqrt(mean), fill=likelihood), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=1/sqrt(true.shape), lty=2) +
    ylab("CV generation") +
    xlab("simulation design") +
    theme(
        legend.position = "none"
    )

coverdf <- df %>%
    group_by(param, likelihood, type) %>%
    summarize(coverage=mean(cover)) %>%
    filter(!(param=="R" & likelihood=="conditional")) %>%
    group_by(param, likelihood, type, coverage) %>%
    mutate(
        lwr=binom.test(100*coverage, 100)$conf.int[1],
        upr=binom.test(100*coverage, 100)$conf.int[2]
    ) %>%
    ungroup %>%
    mutate(
        param=factor(param, levels=c("R", "mean", "shape"), 
                     labels=c("reproductive number", "mean generation", "CV generation"))
    )

g4 <- ggplot(coverdf) +
    geom_point(aes(type, coverage, col=likelihood), position=position_dodge(0.5), size=2) +
    geom_errorbar(aes(type, col=likelihood, ymin=lwr, ymax=upr), position=position_dodge(0.5), width=0.4, lwd=1.1) +
    geom_hline(yintercept=0.95, lty=2) +
    xlab("simulation design") +
    facet_wrap(~param, nrow=3) +
    theme(
        strip.background = element_blank(),
        legend.position=c(0.3, 0.8)
    )

gfinal <- arrangeGrob(
    g1, g2, g3, g4,
    layout_matrix=cbind(c(1, 2, 3), c(4, 4, 4))
)

ggsave("compare.pdf", gfinal, width=6, height=8)
