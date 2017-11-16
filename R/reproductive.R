##' how many people are infected in generation k?
##' @param x simulation result
##' @return vector of number of people in each generation
infection.generation <- function(x) {
    infection <- which(x$t_infected==0)
    generation <- c(length(infection))
    while (TRUE) {
        infection <- which(x$infected_by %in% infection)
        if(length(infection) == 0) {
            break;
        } else {
            generation <- c(generation, length(infection))
        }
    }
    
    generation
}

##' @param x simulation result
##' @param reference reference generation (first generation when there are n infected people)
empirical.R0 <- function(x,
                         reference=100) {
    gen <- infection.generation(x)
    k <- head(which(gen > reference), 1)
    sum(gen[2:(k+1)])/sum(gen[1:k])
}

