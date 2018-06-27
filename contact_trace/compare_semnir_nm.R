library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
library(bbmle)
source("ebola_param_seminr.R")

m.sigma <- 4
n.gamma <- 4

set.seed(101)
nsamp <- 1e5

lat <- rgamma(nsamp, m.sigma, sigma*m.sigma)
inf <- rgamma(nsamp, n.gamma, gamma*n.gamma)

gen <- lat + runif(nsamp, 0, sample(inf, prob=inf, replace=TRUE)) 

true.var <- var(gen)

true.mean <- 1/sigma + (1 + 1/n.gamma)/2 * 1/gamma

true.shape <- true.mean^2/true.var

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

load("seminr_mle_mn.rda")

df <- rlist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="nsample") %>%
    # mutate(nsample=factor(nsample, labels=c(100, 500, 1000))) %>%
    rename(likelihood=method)

g1 <- ggplot(df %>% filter(param=="mean")) +
    geom_boxplot(aes(likelihood, mean, fill=likelihood), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=true.mean, lty=2) +
    geom_hline(yintercept=mle.mean, lty=2, col=2, lwd=1.1)
    ylab("mean generation (days)") +
    xlab("Number of samples") +
    theme(
        legend.position=c(0.7, 0.8)
    )

g2 <- ggplot(df %>% filter(param=="shape")) +
    geom_boxplot(aes(nsample, 1/sqrt(mean), fill=likelihood), width=0.4, alpha=0.5, position=position_dodge(0.5)) +
    geom_hline(yintercept=1/sqrt(true.shape), lty=2) +
    geom_hline(yintercept=1/sqrt(mle.shape), lty=2, col=2, lwd=1.1) +
    ylab("CV generation") +
    xlab("Number of samples") +
    ggtitle("") +
    theme(
        legend.position = "none"
    )

coverdf <- df %>%
    select(param, likelihood, nsample, cover, cover.mle) %>%
    gather(key, value, -param, -likelihood, -nsample) %>%
    mutate(key=factor(key, level=c("cover", "cover.mle"), label=c("true", "MLE"))) %>%
    group_by(param, likelihood, nsample, key) %>%
    mutate(value=ifelse(is.na(value), FALSE, value)) %>%
    summarize(coverage=mean(value)) %>%
    filter(param != "R") %>%
    group_by(param, likelihood, key, coverage, nsample) %>%
    mutate(
        lwr=binom.test(100*coverage, 100)$conf.int[1],
        upr=binom.test(100*coverage, 100)$conf.int[2]
    )

g3 <- ggplot(coverdf) +
    geom_point(aes(nsample, coverage, col=likelihood), position=position_dodge(0.5), size=2) +
    geom_errorbar(aes(nsample, col=likelihood, ymin=lwr, ymax=upr), position=position_dodge(0.5), width=0.2, lwd=1.1) +
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
    g1, g2, g3,
    nrow=1
)

ggsave("compare_seminr_mn.pdf", gfinal, width=8, height=6)
