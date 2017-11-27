library(bbmle)
library(igraph)
library(tidyr)
library(dplyr)

source("../sim/cmp_param.R")
source("../R/generation.R")
source("../R/empirical.R")
source("../R/mle.R")

load("../data/cmp_sim.rda")

addline_format <- function(x,...){
    gsub('\\s','\n',x)
}

cmpGraph <- read.table("../data/CA-CondMat.txt")
cmpGraph <- graph.data.frame(cmpGraph)
cmpGraph <- as.undirected(cmpGraph)
cmpGraph <- delete_vertex_attr(cmpGraph, "name")

d <- degree(cmpGraph)

kappa <- var(d)/mean(d) + mean(d) -1

censor.gi <- lapply(
    reslist
    , network.generation
    , plot=FALSE
    , interval=c(10, 400)
    , interval.type="cases"
)

empirical <- (
    lapply(reslist, empirical.R0, n=100)
    %>% lapply(function(x) data.frame(empirical=x))
    %>% bind_rows(.id="sim")
)

r <- lapply(
    reslist
    , function(x) {
        data <- x$data
        fdata <- data %>% filter(infected <= 400, infected >= 200)
        ll <- lm(log(infected)~time, data=fdata) 
        data.frame(r=ll$coefficients[2])
    }
)

######################################################################

findr <- function(r) (gamma+r)*(sigma+r)/((kappa-1)*sigma-r)-beta

true.r <- uniroot(findr, c(0, 2))$root
mean.r <- mean(unlist(r))

R0fun <- function(r) (1+r/gamma)*(1+r/sigma)
intrinsic.R0 <- R0fun(true.r)
intrinsic.R0_est <- R0fun(mean.r)

R0fun.network <- function(r) (gamma+r)/(gamma*sigma/(sigma+r)+r/kappa)
network.R0 <- R0fun.network(true.r)
network.R0.est <- R0fun.network(mean.r)

local_fun <- function(tau) {
    sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau))
}

localCens_fun <- function(tau) {
    network.R0.est*sigma*(gamma+beta)/(gamma+beta-sigma) * (exp(-sigma*tau)-exp(-(gamma+beta)*tau)) * exp(-mean.r*tau)
}

intrinsic_fun <- function(tau) {
    sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

cens_fun <- function(tau) {
    intrinsic.R0_est*sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau)) * exp(-mean.r*tau)
}

######################################################################

generation <- (
    censor.gi 
    %>% lapply(function(x) data.frame(interval=x)) 
    %>% bind_rows(.id="sim") 
    %>% merge(r %>% bind_rows(.id="sim")) 
    %>% as.tbl 
    %>% group_by(sim)
    %>% mutate(weight=exp(r*interval)/sum(exp(r*interval)))
)

R0.group <- data.frame(
    key=c("contact\ntracing", "temporal\ncorrection", "empirical", "local\ncorrection", "intrinsic"),
    group=factor(c("population based", "population based", "empirical", "individual based", "individual based"),
                 levels=c("population based", "empirical", "individual based"))
)

R0 <- (
    generation 
    %>% group_by(sim) 
    %>% summarize(
        r=mean(r)
        , uncorrected=1/mean(exp(-r*interval))
        , corrected=mean(exp(r*interval)))
    %>% merge(empirical)
    %>% mutate(intrinsic=R0fun(r),
               network=R0fun.network(r))
    %>% gather(key, value, -sim, -r)
    %>% mutate(key=factor(key 
                          , labels=c("contact\ntracing", "temporal\ncorrection", "empirical", "local\ncorrection", "intrinsic")
                          , levels=c("uncorrected", "corrected", "empirical", "network", "intrinsic")))
    %>% as.tbl
    %>% merge(R0.group)
)
empty.df <- data.frame(
    x=c(0.1, 0.1)
    , y=c(0.1, 0.1)
    , group=c("local", "intrinsic")
) %>%
    mutate(group=factor(group, level=.$group))
