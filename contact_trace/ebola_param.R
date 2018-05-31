beta <- 2/5
sigma <- 1/11.4
gamma <- 1/5
N <- 60000

true.mean <- 1/sigma + 1/gamma

intrinsic_fun <- function(tau) {
    sigma*gamma/(gamma-sigma) * (exp(-sigma*tau)-exp(-gamma*tau))
}

true.var <- integrate(function(x) intrinsic_fun(x) * (x-true.mean)^2, lower=0, upper=1000)[[1]]

true.shape <- true.mean^2/true.var

true.R <- beta/gamma