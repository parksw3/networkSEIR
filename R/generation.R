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
                            imin=10,
                            imax=100
                            ) {
    index <- which(!is.na(x$t_infected))
    t_infected <- x$t_infected[index]
    
    tday <- table(floor(t_infected))
    
    day <- as.numeric(names(tday))
    
    case <- as.vector(tday)
    
    firstday <- day[head(which(case >= imin),1)]
    lastday <- day[head(which(case >= imax),1)]
    
    index <- index[t_infected >= firstday & t_infected <= (lastday+1)]
    
    index <- index[order(x$t_infected[index])]
    
    df <- data.frame(
        index=index,
        t_infected=x$t_infected[index],
        case=1:length(index),
        gen=NA
    )
    
    df$gen <- x$t_infected[index] - x$t_infected[x$infected_by[index]]
    
    df$day <- floor(df$t_infected)
    
    df
}
