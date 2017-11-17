library(ggplot2); theme_set(theme_bw())
library(gridExtra)

arrowStep <- function(t, i, astart, adist, asteps, ss){
    pts <- (astart + adist*(0:asteps)/asteps)*ss
    x <- t[pts]
    y <- i[pts]
    xdf <- data.frame(
        x0=x[-length(x)], y0=y[-length(y)],
        x1=x[-1], y1=y[-length(y)],
        type="x"
    )
    
    ydf <- data.frame(
        x0=x[-1], y0=y[-length(y)],
        x1=x[-1], y1=y[-1],
        type="y"
    )
    rbind(xdf, ydf)
}

finTime <- 8
i0 <- 10
C <- 4

steps <- 2
mult <- 6
ss <- mult*steps

df <- data.frame(
    t=(0:(finTime*ss))/ss,
    incidence=i
)

astart <- 3
adist <- 4

segdf <- (
    lapply(
        c(3, 2)
        , arrowStep
        , astart=3
        , adist=4
        , t=t
        , i=i
        , ss=ss) %>%
    bind_rows(.id="type")
)

gg_step <- ggplot(df) +
    geom_line(aes(t, incidence), lwd=1.1) +
    geom_segment(data=segdf, aes(x0, y0, xend=x1, yend=y1), 
                 arrow = arrow(length = unit(0.03, "npc"), type="closed"), lty=2) +
    scale_x_continuous("Time (weeks)", expand=c(0,0), breaks=c(0, 2, 4, 6)) +
    scale_y_continuous("Weekly incidence", expand=c(0,0)) +
    facet_grid(~type) +
    theme(
        panel.grid=element_blank(),
        strip.background = element_blank(),
        strip.text=element_blank()
    )

ggsave("steps.pdf", gg_step, width=8, height=4)
