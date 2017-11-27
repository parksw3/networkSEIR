library(ggplot2); theme_set(theme_bw())
library(gridExtra)

dev.off()

intc <- "#66a61e"

gg1 <- (
    gg_base
    + geom_histogram(
        aes(interval, y=..density..)
        , col='black', fill=obsc
        , alpha=0.5, boundary=0, bins=30) 
    # + geom_line(data=empty.df, aes(x, y, col=group), lwd=1.2)
    # + stat_function(fun=localCens_fun, lwd=1.2, lty=1, xlim=c(0,10), col=localc)
    # + stat_function(fun=cens_fun, lwd=1.2, lty=1, xlim=c(0,10), col=intc)
    + scale_color_manual(values=c(localc, obsc))
    + ggtitle("Observed intervals")
)

gg2 <- (
    gg_base
    + geom_histogram(
        aes(interval, y=..density.., weight=weight)
        , fill=effc, col='black'
        , alpha=0.5, boundary=0, bins=30)
    + geom_line(data=empty.df, aes(x, y, col=group), lwd=1.2)
    + stat_function(fun=local_fun, lwd=1.2, lty=1, xlim=c(0,10), col=localc)
    + stat_function(fun=intrinsic_fun, lwd=1.2, lty=1, xlim=c(0,10), col=intc)
    + ggtitle("Temporal correction")
    + scale_color_manual(values=c(localc, intc))
    + theme(
		 legend.position = c(0.85, 0.85)
		 , legend.title = element_blank()
    )
)

ggsave(pdfname
	, arrangeGrob(gg1, gg2, nrow=1)
	, width=7, height=2.5
)

