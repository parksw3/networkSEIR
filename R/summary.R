summarize.epidemic <- function(data){
    
}



# 
# R0.gen <- list()
# R0.gen[[1]] <- which(infector %in% initial_I)
# j <- 1
# repeat{
#     j <- j + 1
#     R0.gen[[j]] <- which(infector %in% R0.gen[[j-1]])
#     if(length(R0.gen[[j-1]]) > 75){
#         break
#     }
# }
# R0.vec <- unlist(lapply(R0.gen, length))
# R0.generation <- sum(R0.vec[-1])/sum(R0.vec[-length(R0.vec)])
# 
# r.data <- data.frame(
#     tvec <- output[,1],
#     tot <- rowSums(output[,c("I", "R")])
# )
# r.lm <- lm(log(tot)~tvec, data = r.data, tot > 200 & tot < 400)
# r <- unname(coef(r.lm)[2])
# 
# mean.g <- mean(degree(g))
# var.g <- var(degree(g))
# kappa <- var.g/mean.g + mean.g - 1
# 
# R0.Markov <- unname((gamma+r)/(gamma * sigma/(sigma + r) + r/kappa))
# 
# g.summary <- list(R0.generation = R0.generation,
#                 R0.Markov = R0.Markov,
#                 r = r)
# 
# generation <- infected_time - infected_time[infector]
# generation <- generation[!is.na(generation)]
# 
# return(list(result = output,
#             generation = generation,
#             summary = g.summary))