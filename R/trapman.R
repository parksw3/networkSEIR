## Recreate figure 5 from Trapman et al paper
library(igraph)
library(reshape2)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)
source("parameters.R")
windowsFonts(Times=windowsFont("Times New Roman"))
cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

load("../data/condmat_sim.rda")

df <- as.data.frame(do.call("rbind", lapply(sumList, unlist)))
cut <- quantile(df[,1], c(0.05, 0.95))
df2 <- subset(df, little_r > cut[1] & little_r < cut[2])

mL <- melt(df2)

kappa <- var(degree(cmpGraph))/mean(degree(cmpGraph)) + mean(degree(cmpGraph)) - 1

xint <- with(as.list(pars), kappa * beta/(beta + gamma))

theme_custom <- function(){
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(family = "Times"),
        axis.title.y = element_text(family = "Times")
    )
}

g_density <- ggplot(subset(mL, !(variable == "little_r")), aes(value, col = variable, lty = variable)) + 
    geom_line(stat="density") +
    geom_segment(aes(x = xint, y = 0, xend = xint, yend = 2), col = "orange") +
    geom_hline(aes(yintercept = 0), col = "gray") +
    scale_x_continuous(name = expression(basic~reproduction~number~italic(R)[0]),
        breaks = seq(1.4, 2.4, 0.2),
        expand = c(0,0)) +
    scale_colour_manual(values = c("black", "blue", "red")) +
    theme(legend.position = "none") +
    theme_custom()

df3 <- transform(df2,
    ratio = homogeneous_R0/generation_R0)

g_box <- ggplot(df3, aes(x = "", y = ratio)) + 
    stat_boxplot(geom = "errorbar", width = 0.1) +
    geom_boxplot(
        width = 0.3) +
    xlab("") +
    scale_y_continuous(name = expression(ratio~of~estimates~of~italic(R)[0]),
        limit = c(0.75, 1.9),
        breaks = seq(0.8, 1.8, 0.2)
    ) +
    theme(axis.ticks.x = element_blank()) +
    theme_custom()

grid.arrange(g_density, g_box, nrow = 1, widths = c(4,1))
