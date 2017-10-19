##' @param g igraph object
##' @param beta contact rate
##' @param gamma recovery rate
seir.heap <- function(g,
                      beta,
                      gamma,
                      I0=10,
                      seed = NULL){
    if(class(g) != "igraph"){
        stop("g must be an igraph object")
    }
    
    if (!is.null(seed)) set.seed(seed)
    
    V <- as.vector(V(g))
    initial_infected <- sample(V, I0)
    
    t <- 0
    
    queue_v <- queue_t <- queue_infector <- rep(NA, I0)
    
    queue_v[1:I0] <- initial_infected 
    queue_t[1:I0] <- 0
    
    t_infected <- t_recovered <- rep(NA, length(V))
    t_infected[initial_infected] <- 0
    
    t_gillespie <- 0
    c_infected <- I0
    
    intrinsic_generation <- vector('list', length(V))
    
    done <- rep(FALSE, length(V))
    infected_by <- rep(NA, length(V))
    
    target <- initial_infected
    
    i <- I0+1
    
    t <- 0
    
    for (j in initial_infected) {
        n <- as.vector(neighbors(g, j))
        
        rate <- length(n) * beta + gamma
        
        prob <- gamma/rate
        
        ncontact <- rnbinom(1, size=1, prob=prob)
        if (ncontact > 0) {
            queue_v[i:(i+ncontact-1)] <- sample(n, ncontact, replace=TRUE)
            queue_infector[i:(i+ncontact-1)] <- j
        }
        
        time_between <- rexp(ncontact+1, rate=rate)
        intrinsic_generation[[j]] <- c_time <- cumsum(time_between)
        if (ncontact > 0) {
            queue_t[i:(i+ncontact-1)] <- t + c_time[1:ncontact]
        }
        
        t_recovered[j] <- t+tail(c_time,1)
        done[j] <- TRUE
        i <- i+ncontact
    }
    
    stop <- all(done[queue_v])
    
    while (!stop) {
        new_order <- order(queue_t)
        queue_v <- queue_v[new_order]
        queue_t <- queue_t[new_order]
        queue_infector <- queue_infector[new_order]
        
        j.index <- head(which(!done[queue_v]), 1)
        j <- queue_v[j.index]
        infected_by[j] <- queue_infector[j.index]
        t_infected[j] <- queue_t[j.index]
        
        t <- queue_t[j.index]; t_gillespie <- c(t_gillespie, t)
        c_infected <- c_infected +1
        
        n <- as.vector(neighbors(g, j))
        
        rate <- length(n) * beta + gamma
        
        prob <- gamma/rate
        
        ncontact <- rnbinom(1, size=1, prob=prob)
        if (ncontact > 0) {
            queue_v[i:(i+ncontact-1)] <- sample(n, ncontact, replace=TRUE)
            queue_infector[i:(i+ncontact-1)] <- j
        }
        
        time_between <- rexp(ncontact+1, rate=rate)
        intrinsic_generation[[j]] <- c_time <- cumsum(time_between)
        if (ncontact > 0) {
            queue_t[i:(i+ncontact-1)] <- t + c_time[1:ncontact]
        }
        
        t_recovered[j] <- t+tail(c_time,1)
        done[j] <- TRUE
        i <- i+ncontact
        
        stop <- (c_infected == length(V) || all(done[queue_v]))
    }
    
    return(
        list(
            data=data.frame(
                time=t_gillespie,
                infected=I0:c_infected
            ),
            intrinsic_generation=intrinsic_generation,
            t_infected=t_infected,
            t_recovered=t_recovered,
            infected_by=infected_by
        )
    )
}
