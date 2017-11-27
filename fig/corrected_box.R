library(ggplot2); theme_set(theme_bw())
library(gridExtra)

dev.off()

ggsave(pdfname
	, gg_R
	, width=5, height=2.5
)

