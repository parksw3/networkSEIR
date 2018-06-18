library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(bbmle)
source("ebola_param_seminr.R")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

reslist <- vector('list', 2)
names(reslist) <- c("full", "underreport")

load("seminr_mle.rda")

reslist[[1]] <- rlist[[3]]

load("seminr_mle_underreport.rda")

reslist[[2]] <- rlist[[1]]

df <- reslist  %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="type") %>%
    rename(likelihood=method)

g1 <- ggplot(df %>% filter(param=="mean", mean < 1e7)) +
    geom_boxplot(aes(type, mean, fill=likelihood), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.mean, lty=2) +
    ylab("mean generation (days)") +
    xlab("Type") +
    theme(
        legend.position=c(0.4, 0.9)
    )

g2 <- ggplot(df %>% filter(param=="shape")) +
    geom_boxplot(aes(type, 1/sqrt(mean), fill=likelihood), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=1/sqrt(true.shape), lty=2) +
    geom_hline(yintercept=1/sqrt(mle.shape), lty=2, col=2, lwd=1.1) +
    ylab("CV generation") +
    xlab("Type") +
    ggtitle("") +
    theme(
        legend.position = "none"
    )

coverdf <- df %>%
    select(param, likelihood, type, cover, cover.mle) %>%
    gather(key, value, -param, -likelihood, -type) %>%
    mutate(key=factor(key, level=c("cover", "cover.mle"), label=c("true", "MLE"))) %>%
    group_by(param, likelihood, type, key) %>%
    # mutate(value=ifelse(is.na(value), FALSE, value)) %>%
    summarize(coverage=mean(value, na.rm=TRUE),
              n=sum(is.finite(value))) %>%
    # filter(param != "R") %>%
    ungroup %>%
    filter(is.finite(coverage)) %>%
    group_by(param, likelihood, key, coverage, type) %>%
    mutate(
        lwr=binom.test(n*coverage, n)$conf.int[1],
        upr=binom.test(n*coverage, n)$conf.int[2]
    )


g3 <- ggplot(df %>% filter(param=="R", likelihood=="full")) +
    geom_boxplot(aes(type, mean), width=0.4, alpha=0.5, position=position_dodge(0.5), fill="#D95F02") +
    geom_hline(yintercept=true.R, lty=2) +
    ylab("mean generation (days)") +
    xlab("Type") +
    theme(
        legend.position="none"
    )

g4 <- ggplot(coverdf %>% filter(key=="true")) +
    geom_point(aes(type, coverage, col=likelihood), position=position_dodge(0.5), size=2) +
    geom_errorbar(aes(type, col=likelihood, ymin=lwr, ymax=upr), position=position_dodge(0.5), width=0.2, lwd=1.1) +
    geom_hline(yintercept=0.95, lty=2) +
    xlab("Number of samples") +
    scale_y_continuous(limits=c(0,1)) +
    facet_grid(key~param) +
    theme(
        strip.background = element_blank(),
        legend.position = "none",
        panel.spacing = grid::unit(0, "cm")
    )


gfinal <- arrangeGrob(
    g1, g2, g3, g4,
    nrow=2
)

ggsave("compare_seminr_underreport.pdf", gfinal, width=8, height=6)
