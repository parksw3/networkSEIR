library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("ebola_param.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

mlelist <- vector('list', 2)

load("poisson_boot.rda")

mlelist[[1]] <- rlist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="nsample") %>%
    mutate(nsample=factor(nsample, labels=c(100, 500)))

load("seir_boot.rda")

mlelist[[2]] <- rlist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="nsample") %>%
    mutate(nsample=factor(nsample, labels=c(100, 500)))

df <- mlelist %>%
    bind_rows(.id="type") %>%
    mutate(type=factor(type, labels=c("gamma", "SEIR")))

g1 <- ggplot(df %>% filter(param=="mean")) +
    geom_boxplot(aes(nsample, estimate, fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.mean, lty=2) +
    ylab("mean generation (days)") +
    xlab("Number of samples") +
    ggtitle("Gamma") +
    theme(
        legend.position=c(0.7, 0.8)
    )

g2 <- ggplot(df %>% filter(param=="shape")) +
    geom_boxplot(aes(nsample, 1/sqrt(estimate), fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=1/sqrt(true.shape), lty=2) +
    ylab("CV generation") +
    xlab("Number of samples") +
    ggtitle("") +
    theme(
        legend.position = "none"
    )

g3 <- ggplot(df %>% filter(param=="R")) +
    geom_boxplot(aes(nsample, estimate, fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.R, lty=2) +
    ylab("Reproductive number") +
    xlab("Number of samples") +
    ggtitle("") +
    theme(
        legend.position = "none"
    )

coverdf <- df %>%
    group_by(param, type, nsample) %>%
    mutate(cover=ifelse(is.na(cover), FALSE, cover)) %>%
    summarize(coverage=mean(cover)) %>%
    group_by(param, type, nsample, coverage) %>%
    mutate(
        lwr=binom.test(100*coverage, 100)$conf.int[1],
        upr=binom.test(100*coverage, 100)$conf.int[2]
    )

g4 <- ggplot(coverdf) +
    geom_point(aes(nsample, coverage, col=type), position=position_dodge(0.5), size=2) +
    geom_errorbar(aes(nsample, col=type, ymin=lwr, ymax=upr), position=position_dodge(0.5), width=0.2, lwd=1.1) +
    geom_hline(yintercept=0.95, lty=2) +
    xlab("Number of samples") +
    scale_y_continuous(limits=c(0,1)) +
    facet_wrap(~param, nrow=1) +
    theme(
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = grid::unit(0, "cm")
    )

gfinal <- arrangeGrob(
    g1, g2, g3, g4,
    nrow=2
)

ggsave("compare_boot.pdf", gfinal, width=8, height=6)
