##' @param x a vector of nodes
##' @param size number of nodes to pick at random
sample2 <- function(x, size) {
    if(length(x)==1) {
        rep(x, size)
    } else{
        sample(x, size, replace=TRUE)
    }
}

## Simulate an epidemic until done or imax reached
##' @param g igraph object
##' @param beta contact rate
##' @param sigma 1/(latent period)
##' @param gamma 1/(infectious period)
##' @param imax maximum number of infected individuals
seir.heap <- function(g,
                      beta,
                      sigma,
                      gamma,
                      I0,
                      initial_infected,
                      seed = NULL,
                      imax){
    if(class(g) != "igraph"){
        stop("g must be an igraph object")
    }
    
    if (!is.null(seed)) set.seed(seed)
    
    V <- as.vector(V(g))
    
    if (missing(initial_infected)) {
        if (missing(I0)) stop("specify the initial conditions")
        
        initial_infected <- sample(V, I0)    
    }

    I0 <- length(initial_infected)
    
    if (missing(imax)) imax <- length(V)
    
    t <- 0
    
    queue_v <- queue_t <- queue_infector <- rep(NA, I0)
    
    queue_v[1:I0] <- initial_infected 
    queue_t[1:I0] <- 0
    
    t_infected <- t_infectious <- t_recovered <- rep(NA, length(V))
    t_infected[initial_infected] <- 0
    
    t_gillespie <- 0
    c_infected <- I0
    
    intrinsic_generation <- vector('list', length(V))
    
    done <- rep(FALSE, length(V))
    infected_by <- rep(NA, length(V))
    
    target <- initial_infected
    
    t <- 0
    
    for (j in initial_infected) {
        latent <- rexp(1, sigma)
        t_infectious[j] <- t + latent
        
        n <- as.vector(neighbors(g, j))
        
        rate <- length(n) * beta + gamma
        
        prob <- gamma/rate
        
        ncontact <- rnbinom(1, size=1, prob=prob)
        if (ncontact > 0) {
            queue_v <- c(queue_v, sample2(n, ncontact))
            queue_infector <- c(queue_infector, rep(j, ncontact))
        }
        
        time_between <- rexp(ncontact+1, rate=rate)
        c_time <- latent + cumsum(time_between)
        intrinsic_generation[[j]] <-  c_time[1:ncontact]
        if (ncontact > 0) {
            queue_t <- c(queue_t, t + c_time[1:ncontact])
        }
        
        t_recovered[j] <- t+tail(c_time,1)
        done[j] <- TRUE
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
        
        latent <- rexp(1, sigma)
        t_infectious[j] <- t + latent
        
        c_infected <- c_infected +1
        
        n <- as.vector(neighbors(g, j))
        
        rate <- length(n) * beta + gamma
        
        prob <- gamma/rate
        
        ncontact <- rnbinom(1, size=1, prob=prob)
        if (ncontact > 0) {
            queue_v <- c(queue_v, sample2(n, ncontact))
            queue_infector <- c(queue_infector, rep(j, ncontact))
        }
        
        time_between <- rexp(ncontact+1, rate=rate)
        c_time <- latent + cumsum(time_between)
        intrinsic_generation[[j]] <- c_time[1:ncontact]
        if (ncontact > 0) {
            queue_t <- c(queue_t, t + c_time[1:ncontact])
        }
        
        t_recovered[j] <- t+tail(c_time,1)
        done[j] <- TRUE
        
        stop <- (c_infected == length(V) || all(done[queue_v]) || c_infected > imax)
    }
    
    return(
        list(
            data=data.frame(
                time=t_gillespie,
                infected=I0:c_infected
            ),
            intrinsic_generation=intrinsic_generation,
            t_infected=t_infected,
            t_infectious=t_infectious,
            t_recovered=t_recovered,
            infected_by=infected_by
        )
    )
}
