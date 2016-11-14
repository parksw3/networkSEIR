seir.gillespie <- function(g, parameters,
                           R0_threshold = 75,
                           r_interval = c(200, 400),
                           verbose = FALSE,
                           seed = NULL){
    if(class(g) != "igraph"){
        stop("g must be an igraph object")
    }
    
    with(as.list(parameters),{
        N <- length(V(g))
        
        state <- rep("S", N)
        rate <- rep(0, N)
        
        if(!is.null(seed)) set.seed(seed)
        initial_I <- sample(N, I0)
        state[initial_I] <- "I"
            
        inf.partner <- rep(0, N)
        
        for(i in initial_I){
            inf.partner[neighbors(g, i)] <- inf.partner[neighbors(g,i)] + 1
        }
        
        rate <- beta * inf.partner
        
        rate[initial_I] <- gamma
        
        infected_time <- rep(NA, N)
        infected_time[initial_I] <- 0
        infector <- rep(NA, N)
        infected_order <- rep(NA, N)
        infected_order[1:I0] <- initial_I
        
        rsum <- sum(rate)
        
        t <- 0; S <- N-I0; E <- 0; I <- I0; R <- 0
        
        output <- vector("list", 20000)
        
        i <- 1
        
        output[[i]] <- data.frame(t, S, E, I, R)
        
        while(rsum > 0){
            i <- i + 1
            next.t <- rexp(1, rsum)
            t <- t + next.t
            
            vchange <- which(rmultinom(1, size = 1, prob = rate) == 1)
            
            if(state[vchange] == "S"){
                state[vchange] <- "E"
                rate[vchange] <- sigma
                
                n.V <- neighbors(g,vchange)
                n.V.I <- n.V[state[n.V] == "I"]
                infector[vchange] <- n.V.I[sample(length(n.V.I), 1)]
                
                S <- S - 1; E <- E + 1
                
                infected_order[N - S] <- vchange
            }else if(state[vchange] == "E"){
                state[vchange] <- "I"
                
                n.V <- neighbors(g,vchange)
                n.V.S <- n.V[state[n.V] == "S"]
                inf.partner[n.V.S] <- inf.partner[n.V.S] + 1
                rate[n.V.S] <- beta * inf.partner[n.V.S]
                rate[vchange] <- gamma
                
                infected_time[vchange] <- t
                
                E <- E - 1; I <- I + 1
            }else if(state[vchange] == "I"){
                state[vchange] <- "R"
                n.V <- neighbors(g, vchange)
                n.V.S <- n.V[state[n.V] == "S"]
                inf.partner[n.V.S] <- inf.partner[n.V.S] - 1
                rate[n.V.S] <- beta * inf.partner[n.V.S]
                rate[vchange] <- 0
                
                I <- I - 1; R <- R + 1
            }
            
            rsum <- sum(rate)
            output[[i]] <- data.frame(t, S, E, I, R)
            if(verbose & i %% 1000 == 0){
                print(output[[i]])
            }
        }
        
        result <- rbindlist(output)
        class(result) <- "data.frame"
        
        data <- list(
            infected_time = infected_time,
            infector = infector,
            infected_order = infected_order
        )
        
        ## calculate exponential growth rate
        r_data <- data.frame(
            tvec = result[,1],
            tot = rowSums(result[,c("I", "R")])
        )
        
        if(max(r_data[,"tot"]) < r_interval[2]){
            stop("epidemic size is too small")
        }else{
            r.lm <- lm(log(tot)~tvec, data = r_data, tot > r_interval[1] & tot < r_interval[2])
            r <- unname(coef(r.lm)[2])
            
            ## calculate R0 based on generation
            genList <- list()
            genList[[1]] <- which(infector %in% initial_I)
            j <- 1
            repeat{
                j <- j + 1
                genList[[j]] <- which(infector %in% genList[[j-1]])
                if(length(genList[[j-1]]) >= R0_threshold){
                    break
                }
            }
            generation_vec <- unlist(lapply(genList, length))
            generation_R0 <- sum(generation_vec[-1])/sum(generation_vec[-length(generation_vec)])
            
            ## calculate R0 based on the exponential growth rate 
            ## under configuration model assumption
            mean.g <- mean(degree(g))
            var.g <- var(degree(g))
            kappa <- var.g/mean.g + mean.g - 1
            network_R0 <- (gamma+r)/(gamma * sigma/(sigma + r) + r/kappa)
            
            ## calculate R0 based on the exponential growth rate 
            ## under homogeneous model assumption
            homogeneous_R0 <- (1 + r/gamma) * (1 + r/sigma)
            
            ## store result in a list
            summary <- list(
                little_r = r,
                generation_R0 = generation_R0,
                network_R0 = network_R0,
                homogeneous_R0 = homogeneous_R0
            )
            
            return(
                list(
                    epidemic.summary = summary,
                    edpidemic.result = result,
                    epidemic.data = data
                )
            )
        }
    })
}
