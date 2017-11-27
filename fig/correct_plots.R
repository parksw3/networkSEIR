library(ggplot2); theme_set(theme_bw())
library(gridExtra)

obsc <- "#e7298a"
effc <- "#7570b3"
localc <- "#d95f02"

gg_base <- (
    ggplot(generation)
    + scale_x_continuous(name="Generation interval") 
    + theme(
        panel.grid = element_blank()
        , panel.border=element_blank()
        , axis.line.x=element_line()
        , axis.title.y = element_blank()
        , axis.text.y = element_blank()
        , axis.ticks.y = element_blank()
        , plot.title = element_text(hjust=0.5)
    )
)

gg_R <- (
    ggplot(R0, aes(key, value, fill=key)) 
    + geom_boxplot(alpha=0.5, width=0.7)
    + scale_y_log10("Reproductive number", breaks=c(2, 5, 10, 20))
    + scale_fill_manual(values=c(obsc, effc, 1, "#d95f02", "#66a61e") )
    + facet_grid(~group, scale="free", space="free_x")
    + theme(
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank(),
        legend.position = "none",
        panel.grid=element_blank(),
        axis.title.x=element_blank()
    )
)

# gg_correction <- arrangeGrob(gg1, gg2, gg_R, nrow=1)

