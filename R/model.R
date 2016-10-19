seir.gillespie <- function(g, parameters,
                           verbose = FALSE){
    if(class(g) != "igraph"){
        stop("g must be an igraph object")
    }
    
    with(as.list(parameters),{
        N <- length(V(g))
        
        state <- rep("S", N)
        rate <- rep(0, N)
        
        infected <- sample(N, I0)
        state[infected] <- "I"
        
        inf.partner <- rep(0, N)
        
        for(i in infected){
            inf.partner[neighbors(g, i)] <- inf.partner[neighbors(g,i)] + 1
        }
        
        rate <- beta * inf.partner
        
        rate[infected] <- gamma
        
        inf.t <- rep(NA, N)
        inf.t[infected] <- 0
        infector <- rep(NA, N)
        inf.order <- rep(NA, N)
        inf.order[1:I0] <- infected
        
        
        rsum <- sum(rate)
        
        t <- 0; S <- N-I0; E <- 0; I <- I0; R <- 0
        
        output <- list()
        
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
                
                inf.order[N - S] <- vchange
            }else if(state[vchange] == "E"){
                state[vchange] <- "I"
                
                n.V <- neighbors(g,vchange)
                n.V.S <- n.V[state[n.V] == "S"]
                inf.partner[n.V.S] <- inf.partner[n.V.S] + 1
                rate[n.V.S] <- beta * inf.partner[n.V.S]
                rate[vchange] <- gamma
                
                inf.t[vchange] <- t
                
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
        
        output <- do.call('rbind', output)
        
        R0.gen <- list()
        R0.gen[[1]] <- which(infector %in% infected)
        j <- 1
        repeat{
            j <- j + 1
            R0.gen[[j]] <- which(infector %in% R0.gen[[j-1]])
            if(length(R0.gen[[j-1]]) > 75){
                break
            }
        }
        R0.vec <- unlist(lapply(R0.gen, length))
        R0.generation <- sum(R0.vec[-1])/sum(R0.vec[-length(R0.vec)])
        
        r.data <- data.frame(
            tvec <- output[,1],
            tot <- rowSums(output[,c("I", "R")])
        )
        r.lm <- lm(log(tot)~tvec, data = r.data, tot > 200 & tot < 400)
        r <- unname(coef(r.lm)[2])
        
        mean.g <- mean(degree(g))
        var.g <- var(degree(g))
        kappa <- var.g/mean.g + mean.g - 1
        
        R0.Markov <- unname((gamma+r)/(gamma * sigma/(sigma + r) + r/kappa))
        
        g.summary <- list(R0.generation = R0.generation,
                        R0.Markov = R0.Markov,
                        r = r)
        
        generation <- inf.t - inf.t[infector]
        generation <- generation[!is.na(generation)]
        
        return(list(result = output,
                    generation = generation,
                    summary = g.summary))
    })
}