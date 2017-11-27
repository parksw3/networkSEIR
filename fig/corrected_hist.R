library(ggplot2); theme_set(theme_bw())
library(gridExtra)

dev.off()

ggsave(pdfname
	, arrangeGrob(gg1, gg2, nrow=1)
	, width=7, height=2.5
)

