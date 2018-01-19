intrinsic.generation <- function(x, plot=TRUE, ...) {
    ig <- unlist(x$intrinsic_generation)
    if(plot) hist(ig, freq=FALSE, ...)
    invisible(ig)
}

network.generation <- function(x, plot=TRUE,
                               type=c("backward", "forward"),
                               interval,
                               interval.type=c("time", "cases"),
                               ...) {
    type <- match.arg(type)
    interval.type <- match.arg(interval.type)
    
    if (missing(interval)) {
        interval <- c(0, max(x$data$time))
    } else if (interval.type=="cases") {
        interval <- x$data$time[match(interval, x$data$infected)]
    }
    
    cohort <- which(x$t_infected >= interval[1] & x$t_infected <= interval[2])
    
    if (type == "backward") {
        infector <- x$infected_by[cohort]
        ng <- x$t_infected[cohort] -  x$t_infected[infector]
        ng <- ng[!is.na(ng)]
    } else {
        ng <- unlist(lapply(cohort, function(z) {
            infectee <- which(z==x$infected_by)
            x$t_infected[infectee] - x$t_infected[z]
        }))
    }
    
    if (plot) hist(ng, freq=FALSE, ...)
    
    invisible(ng)
}

generation.data <- function(x,
                            incidence.report=1,
                            generation.report=1,
                            cutoff=200
                            ) {
    index <- which(!is.na(x$t_infected))
    index <- index[tail(order(x$t_infected[index]), -cutoff)]
    
    incidence.index <- index[runif(length(index)) < incidence.report]
    
    generation.index <- incidence.index[runif(length(incidence.index)) < generation.report]
    
    df <- data.frame(
        index=incidence.index,
        t_infected=x$t_infected[incidence.index],
        case=1:length(incidence.index)+cutoff,
        gen=NA
    )
    
    df$gen[match(generation.index, incidence.index)] <- x$t_infected[generation.index] - x$t_infected[x$infected_by[generation.index]]
    
    df$week <- floor(df$t_infected/7)
    
    df[df$week > min(df$week) & df$week < max(df$week),]
}
