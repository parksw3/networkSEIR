library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("ebola_param_seminr.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

mlelist <- vector('list', 2)

load("kernel_ma.rda")

mlelist[[1]] <- rlist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="nsample")

load("seminr_ma.rda")

mlelist[[2]] <- rlist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="nsample")

df <- mlelist %>%
    bind_rows(.id="type") %>%
    mutate(type=factor(type, labels=c("gamma", "SEmInR")))

g1 <- ggplot(df %>% filter(param=="mean")) +
    geom_boxplot(aes(type, mean, fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.mean, lty=2) +
    ylab("mean generation (days)") +
    xlab("Type") +
    theme(
        legend.position=c(0.2, 0.8)
    )

g2 <- ggplot(df %>% filter(param=="shape")) +
    geom_boxplot(aes(type, 1/sqrt(mean), fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=1/sqrt(true.shape), lty=2) +
    ylab("CV generation") +
    xlab("Type") +
    theme(
        legend.position = "none"
    )

g3 <- ggplot(df %>% filter(param=="r")) +
    geom_boxplot(aes(type, mean, fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    # geom_hline(yintercept=1/sqrt(true.shape), lty=2) +
    ylab("Growth rate") +
    xlab("Type") +
    theme(
        legend.position = "none"
    )

g4 <- ggplot(df %>% filter(param=="R")) +
    geom_boxplot(aes(type, mean, fill=type), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.R, lty=2) +
    ylab("Reproductive number") +
    xlab("Type") +
    theme(
        legend.position = "none"
    )

coverdf <- df %>%
    group_by(param, type) %>%
    filter(!is.na(cover)) %>%
    summarize(coverage=mean(cover)) %>%
    group_by(param, type) %>%
    mutate(
        lwr=binom.test(100*coverage, 100)$conf.int[1],
        upr=binom.test(100*coverage, 100)$conf.int[2]
    ) 

g5 <- ggplot(coverdf) +
    geom_point(aes(type, coverage, col=type), position=position_dodge(0.5), size=2) +
    geom_errorbar(aes(type, col=type, ymin=lwr, ymax=upr), position=position_dodge(0.5), width=0.2, lwd=1.1) +
    geom_hline(yintercept=0.95, lty=2) +
    xlab("Type") +
    scale_y_continuous(limits=c(0,1)) +
    facet_wrap(~param, nrow=2) +
    theme(
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = grid::unit(0, "cm")
    )

gfinal <- arrangeGrob(
    g1, g2, g4, g5
)

ggsave("compare_ma.pdf", gfinal, width=8, height=6)
