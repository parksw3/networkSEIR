library(bbmle)
library(epigrowthfit)
source("ebola_param_seminr.R")
source("../R/generation.R")
source("../R/empirical.R")
source("minuslogl.R")

load("kernel_sim_seminr.rda")

## approximate...
approx.r <- 1/2*(-(sigma+gamma)+sqrt((sigma-gamma)^2 + 4 * beta * sigma))

sens.mean <- function(log.shape, log.mean, log.r) {
    
    unname(c(
        -exp(log.shape+2*log.mean+log.r)/(exp(log.shape) - exp(log.mean+log.r))^2,
        exp(2*log.shape - log.mean)/(exp(log.shape - log.mean) - exp(log.r))^2,
        exp(log.shape+2*log.mean+log.r)/(exp(log.shape) - exp(log.mean+log.r))^2
    ))
}

sens.R <- function(log.shape, log.mean, log.r) {
    
    unname(c(
        (1-exp(-log.shape+log.mean+log.r))^(-exp(log.shape)) * 
            (-exp(log.mean + log.r)/(1-exp(-log.shape+log.mean+log.r))-exp(log.shape) * log(1 - exp(-log.shape+log.mean+log.r)))  ,
        exp(log.mean+log.r) * (1 - exp(-log.shape+log.mean+log.r))^(-exp(log.shape)-1),
        exp(log.mean+log.r) * (1 - exp(-log.shape+log.mean+log.r))^(-exp(log.shape)-1)
    ))
}

rlist <- vector('list', 1)

report <- 0.7

nsample <- c(999)

set.seed(101)
for (j in 1:length(nsample)) {
    
    subrlist <- vector('list', 100)
    
    for (i in 1:100) {
        print(i)
        sim <- reslist[[i]]
        
        tmax <- max(sim$data$time)
        
        data <- generation.data(sim, tmax=tmax, "backward")[1:nsample[j],]
        
        data <- data[runif(nsample[j]) <= report,]
        
        itime <- data$t_infected[!duplicated(data$index)]
        
        incdf <- as.data.frame(table(round(itime)))
        
        incdf$Var1 <- as.numeric(as.character(incdf$Var1))
        
        tvec <- seq(min(incdf$Var1), max(incdf$Var1), by=1)
        
        ivec <- rep(0, length(tvec))
        
        oo <- match(tvec, incdf$Var1)
        
        ivec[!is.na(oo)] <- incdf$Freq[oo[!is.na(oo)]]
        
        ee <- epigrowthfit(
            time=tvec,
            deaths=ivec,
            theta0=c(r=approx.r, x0=exp(-10), K=exp(10), ll.k=20),
            model=get_model("logistic")
        )
        
        if (coef(ee)[4] > 1e10 || coef(ee)[1] < 1e-5) {
            ee <- epigrowthfit(
                time=tvec,
                deaths=ivec,
                theta0=c(r=approx.r, x0=exp(-10), K=exp(7)),
                model=get_model("logistic"),
                distrib="poisson"
            )
        }
        
        gen <- data$t_infected - data$t_infected[match(data$infected_by, data$index)]
        gen <- gen[!is.na(gen)]
        
        m <- mle2(gen~dgamma(shape=exp(log.shape), rate=exp(log.shape)/exp(log.mean)),
             start=list(log.shape=log(true.shape), log.mean=log(true.mean)),
             data=data.frame(gen=gen)
        )
        
        est.r <- coef(ee)[1]
        est.shape <- exp(coef(m)[1])
        
        est.mean <- unname(est.shape/(est.shape/exp(coef(m)[2]) - est.r))
        
        est.R <- (1 - exp(coef(m)[2])/est.shape * est.r)^(-est.shape)
        
        ci.r <- sort(exp(confint(ee@profile)))
        ci <- confint(m)
        ci.shape <- exp(ci[1,])
        
        est.vcov <- matrix(0, 3, 3)
        
        est.vcov[1:2, 1:2] <- vcov(m)
        
        if (vcov(ee@mle2)[1,1]==0) {
            ## extremely approximate....
            est.vcov[3,3] <- (diff(log(ci.r))/(2*1.96))^2
        } else {
            est.vcov[3,3] <- vcov(ee@mle2)[1,1]
        }
        
        var.mean <- unname(c(sens.mean(coef(m)[1], coef(m)[2], coef(ee@mle2)[1]) %*% est.vcov %*% sens.mean(coef(m)[1], coef(m)[2], coef(ee@mle2)[1])))
        
        ci.mean <- c(est.mean-sqrt(var.mean)*1.96, est.mean+sqrt(var.mean)*1.96)
        
        var.R <- unname(c(sens.R(coef(m)[1], coef(m)[2], coef(ee@mle2)[1]) %*% est.vcov %*% sens.R(coef(m)[1], coef(m)[2], coef(ee@mle2)[1])))
        
        ci.R <- c(est.R-sqrt(var.R)*1.96, est.R+sqrt(var.R)*1.96)
        
        res1 <- data.frame(
            method="Ma",
            param=c("r", "R", "mean", "shape"),
            mean=c(est.r, est.R, est.mean, est.shape),
            upr=c(ci.r[2], ci.R[2], ci.mean[2], ci.shape[2]),
            lwr=c(ci.r[1], ci.R[1], ci.mean[1], ci.shape[1]),
            cover=c(
                NA,
                ci.R[1] < true.R && true.R < ci.R[2],
                ci.mean[1] < true.mean && true.mean < ci.mean[2],
                ci.shape[1] < true.shape && true.shape < ci.shape[2]
            )
        )
        
        rownames(res1) <- NULL
        
        print(res1)
        
        subrlist[[i]] <- res1
    }
    
    rlist[[j]] <- subrlist
}    

save("rlist", file="kernel_ma.rda")
