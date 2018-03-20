library(bbmle)
library(tidyr)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("../sim/full_param2.R")
source("../R/generation.R")

load("R0_sim.rda")

t.interval <- seq(50, 350, by=10)

estlist <- vector('list', 10)
simlist <- vector('list', 10)

for (i in 1:10) {
    res <- reslist[[i]]
    
    data <- generation.data(res, imin=10, imax=100)
    day <- sort(unique(data$day))
    case <- as.vector(table(data$day))
    
    r.start <- unname(coef(lm(log(case)~day))[2])
    
    fit <- bbmle::mle2(
        case~dnbinom(mu=exp(log.i0)*exp(r*day), size=exp(log.size))
        , start=list(r=r.start, log.i0=0, log.size=0)
        , data=data.frame(
            case=case
            , day=sort(unique(data$day))
        )
        , method="BFGS"
    )
    
    r <- coef(fit)[["r"]]
    
    estimate <- lapply(t.interval, function(t) {
        gg <- network.generation(reslist[[i]], interval=c(0,t), plot=FALSE)
        data.frame(
            time=t,
            generation=mean(gg),
            R0=1/mean(exp(-r*gg))
        )
    })
    
    estlist[[i]] <- do.call("rbind", estimate)
    
    simlist[[i]] <- data.frame(
        time=sort(unique(round(res$data$time))),
        cases=as.vector(table(round(res$data$time)))
    )
}

simdf <- simlist %>%
    bind_rows(.id="sim") %>%
    group_by(time) %>%
    summarize(
        mean=mean(cases),
        lwr=quantile(cases, 0.025),
        upr=quantile(cases, 0.975)
    )

gsim <- ggplot(simdf, aes(time, mean)) +
    geom_ribbon(aes(ymin=lwr, ymax=upr), alpha=0.1) +
    geom_line() +
    ggtitle("Daily incidence") +
    scale_x_continuous("Days", limits=c(0, 350)) +
    scale_y_continuous("Cases")
    
estdf <- estlist %>%
    bind_rows(.id="sim") %>%
    gather(key, value, -sim, -time) %>%
    mutate(key=factor(key, labels=c("Mean generation interval (days)", "Reproductive number")) )

gest <- ggplot(estdf) +
    geom_boxplot(aes(time,value, group=time)) +
    facet_wrap(~key, scale="free_y", ncol=1) +
    scale_x_continuous("Days", limits=c(0, 350)) +
    theme(
        strip.background = element_blank(),
        axis.title.y = element_blank()
    )

g1 <- arrangeGrob(gsim, gest, ncol=2)

ggsave("figure1.pdf", g1, width=10, height=6)

