library(ggplot2); theme_set(theme_bw())
library(gridExtra)

dev.off()

ggsave(pdfname
	, gg_R
	, width=4, height=2.5
)

