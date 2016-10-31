generation <- function(data,
                       interval = c(0, 10),
                       type = "forward"){
    with(data,{
        cohort <- which(infected_time > interval[1] & infected_time < interval[2])
        if(type == "forward"){
            cohort.i <- which(infector %in% cohort)
            gen <- infected_time[cohort.i] - infected_time[infector[cohort.i]]
        }else if(type == "backward"){
            cohort.i <- infector[cohort]
            gen <- infected_time[cohort] - infected_time[cohort.i]
        }
        return(gen)
    })
}