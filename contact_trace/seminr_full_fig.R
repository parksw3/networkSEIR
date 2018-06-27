library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
source("ebola_param_seminr.R")

load("seminr_full_sim.rda")
load("seminr_full_mle.rda")

incdf <- reslist %>%
    lapply("[[", "data") %>%
    lapply(function(x){
        tb <- table(round(x$time))
        
        data.frame(
            day=as.numeric(names(tb)),
            cases=as.numeric(tb)
        )
    }) %>%
    bind_rows(.id="sim")

mledf <- rlist %>%
    lapply(bind_rows, .id="sim") %>%
    bind_rows(.id="time") %>%
    mutate(time=factor(time, labels=c(50, 100, 150, 200)))

meandf <- mledf %>%
    filter(param=="mean") %>%
    filter(mean <1e5)


ggplot(meandf) +
    geom_boxplot(aes(time, mean, col=method)) +
    geom_hline(yintercept=true.mean, lty=2)





