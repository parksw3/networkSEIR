source("../sim/full_param.R")
library(deSolve)
library(dplyr)
library(ggplot2); theme_set(theme_bw())
library(gridExtra)

R0fun <- function(r) (1+r/gamma)*(1+r/sigma)

R0fun.network <- function(r, kappa) (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)

genfun.network <- function(r, kappa) {
    if (is.na(kappa)) {
        1/sigma + 1/gamma
    } else {
        beta <- (gamma+r) * (sigma + r)/((kappa-1)*sigma-r)
        1/sigma + 1/(gamma+beta)
    }
}

r <- seq(0, 4, 0.01)

kappa <- c(5, 10, 50, 100, 500)

h.R0 <- data.frame(
    r=r,
    R0=R0fun(r),
    kappa=NA,
    type="homogeneous"
)

n.R0 <- lapply(kappa,
    function(k) data.frame(
        r=r,
        R0=R0fun.network(r,k),
        type="network"
    )
)

names(n.R0) <- kappa

generation <- n.R0 %>%
    bind_rows(.id="kappa") %>%
    rbind(h.R0) %>%
    as.tbl %>%
    mutate(kappa=as.numeric(kappa)) %>%
    mutate(gi=genfun.network(r, kappa)) %>%
    mutate(rho=r*gi)

ggplot(generation, aes(r, R0, group=kappa)) +
    geom_line()
