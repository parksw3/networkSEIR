estfun <- function(data) {
    unrecorded <- !(data$infected_by %in% data$index)
    
    unobs <- data.frame(
        index=data$infected_by[unrecorded],
        infected_by=NA,
        t_infected=data$t_infected[unrecorded] - data$generation[unrecorded],
        generation=NA
    )
    unobs <- unobs[!duplicated(unobs),]
    
    data <- rbind(data, unobs)
    
    tmax <- max(data$t_infected)
    
    data$weight <- sapply(data$generation, function(x) 1/mean(tmax - data$t_infected >= x)) 
    
    ## I think this is the natural estimator for R?
    R <- mean(data$weight, na.rm=TRUE)
    
    mean.gen <- weighted.mean(data$generation, w=data$weight, na.rm=TRUE)
    
    var.gen <- weighted.mean((data$generation-mean.gen)^2, w=data$weight, na.rm=TRUE)
    
    data.frame(
        mean=mean.gen,
        shape=mean.gen^2/var.gen,
        R=R
    )
}

bootfun <- function(data, nrep=500) {
    bootlist <- replicate(nrep, estfun(data[sample(1:nrow(data), replace=TRUE),]), simplify=FALSE)

    bootdata <- do.call("rbind", bootlist)
    
    cidata <- apply(bootdata, 2, quantile, c(0.025, 0.975))
    
    est <- estfun(data)

    data.frame(
        method="bootstrap",
        param=c("mean", "shape", "R"),
        estimate=unlist(est),
        lwr=cidata[1,],
        upr=cidata[2,]
    )
}
