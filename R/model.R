seir.gillespie <- function(g, parameters,
                           verbose = FALSE){
    if(class(g) != "igraph"){
        stop("g must be an igraph object")
    }
    
    with(as.list(parameters),{
        N <- length(V(g))
        state <- rep("S", N)
        
        connected <- as(g[], "sparseMatrix")
        connected <- unname(connected)
        ## connected <- as.matrix(g[])
        max.edge <- which.max(colSums(connected))
        state[max.edge] <- "I"
        
        inf.t <- rep(NA, N)
        inf.t[max.edge] <- 0
        infector <- rep(NA, N)
        inf.order <- rep(NA, 0)
        
        inf <- as.numeric(state == "I")
        
        rate.mat <- beta * inf * connected
        rate.mat[max.edge,max.edge] <- gamma
        
        rsum <- sum(rate.mat)
        
        t <- 0; S <- N-1; E <- 0; I <- 1; R <- 0
        
        output <- list()
        
        i <- 1
        
        output[[i]] <- data.frame(t, S, E, I, R)
        
        while(rsum > 0){
            i <- i + 1
            next.t <- rexp(1, rsum)
            t <- t + next.t
            
            # prob <- c(rate.mat)
            
            prob.ind <- which(rate.mat != 0, arr.ind = TRUE)
            prob <- rate.mat[prob.ind]
            
            next.event <- which(rmultinom(1, size = 1, prob = prob) == 1)
            change.ind <- prob.ind[next.event,]
            vchange <- change.ind[2]
            
            if(state[vchange] == "S"){
                state[vchange] <- "E"
                rate.mat[,vchange] <- 0
                rate.mat[vchange,vchange] <- sigma
                infector[vchange] <- change.ind[1]
                    
                S <- S - 1; E <- E + 1
                
                inf.order[N - S] <- vchange
            }else if(state[vchange] == "E"){
                state[vchange] <- "I"
                rate.mat[vchange,(state == "S")] <- beta
                rate.mat[vchange,vchange] <- gamma
                inf.t[vchange] <- t
                
                E <- E - 1; I <- I + 1
            }else if(state[vchange] == "I"){
                state[vchange] <- "R"
                rate.mat[vchange,] <- 0
                
                I <- I - 1; R <- R + 1
            }
            
            rsum <- sum(rate.mat)
            output[[i]] <- data.frame(t, S, E, I, R)
            if(verbose){
                print(output[[i]])
            }
        }
        
        output <- do.call('rbind', output)
        
        generation <- inf.t - inf.t[infector]
        generation <- generation[!is.na(generation)]
        
        return(list(result = output,
                    generation = generation))  
    })
}
