library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("ebola_param_seminr.R")

load("seminr_full_sim.rda")
load("seminr_full_mle_example.rda")

scale_colour_discrete <- function(...,palette="Dark2") scale_colour_brewer(...,palette=palette)
scale_fill_discrete <- function(...,palette="Dark2") scale_fill_brewer(...,palette=palette)

mledf <- rlist %>%
    bind_rows(.id="time") %>%
    mutate(time=factor(time, levels=1:37, labels=seq(20, 200, by=5))) %>%
    mutate(time=as.numeric(as.character(time)))

meandf <- mledf %>%
    filter(param=="mean")

gmean <- ggplot(meandf) +
    geom_line(aes(time, mean, col=method)) +
    geom_ribbon(aes(time, ymin=lwr, ymax=upr, fill=method, col=method), alpha=0.1) +
    scale_x_continuous(limits=c(0, 250), expand=c(0,0)) +
    ylab("GI mean") +
    coord_cartesian(ylim = c(8, 25)) +
    geom_hline(yintercept=true.mean, lty=2) +
    theme(
        legend.position="none",
        legend.title=element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 15, b = 0, l = 0))
    )

shapedf <- mledf %>%
    filter(param=="shape")

gshape <- ggplot(shapedf) +
    geom_line(aes(time, mean, col=method)) +
    geom_ribbon(aes(time, ymin=lwr, ymax=upr, fill=method, col=method), alpha=0.1) +
    geom_hline(yintercept=true.shape, lty=2) +
    ylab("GI shape") +
    scale_x_continuous(limits=c(0, 250), expand=c(0,0), breaks=c(0, 50, 100, 150, 200)) +
    theme(
        legend.position=c(0.8, 0.8),
        legend.title=element_blank(),
        legend.direction = "horizontal"
    )

df <- reslist[[1]]$data
tb <- table(round(df$time))
incdf <- data.frame(
    day=as.numeric(names(tb)),
    cases=as.vector(tb)
)

ginc <- ggplot(incdf) +
    geom_line(aes(day, cases)) +
    scale_x_continuous(limits=c(0, 250), expand=c(0,0)) +
    scale_y_log10("daily incidence") +
    theme(
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_text(margin = margin(t = 0, r = 9, b = 0, l = 0))
    )

gex <- arrangeGrob(
    ginc, gmean, gshape,
    nrow=3
)

ggsave("example.pdf", gex, width=8, height=6)
