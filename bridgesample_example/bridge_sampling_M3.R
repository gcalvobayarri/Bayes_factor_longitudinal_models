# devtools::install_github("quentingronau/bridgesampling@master")

# install.packages("bridgesampling")
# install.packages("BayesFactor")
# install.packages("rjags")
# install.packages("R2jags")
# install.packages("rstan")
# install.packages("sm")

library("bridgesampling")
library("BayesFactor")
library("R2jags")
library("rstan")
library("sm")

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ Controversial data @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#-------------------------------------------------------------------------------
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

load("./bridgesample_example/long_data_controversial.RData")
y <- Y
rm(Y)
n <- length(y[, 1])
T <- length(y[1, ])

sizes <- c()
for(i in 1 : n){
  sizes[i] <- T - sum(is.na(y[i,]))
}

#-------------------------------------------------------------------------------
# JAGS models
#-------------------------------------------------------------------------------

code_H3 <-
  "model{
  for (i in 1 : n){
    y[i, 1] ~ dnorm(mu[i, 1], lambda0)
    mu[i, 1] <- beta0 + b1[i] * times[i, 1]
    for (t in 2 : sizes[i]){
      y[i, t] ~ dnorm(Mi[i, t], lambda) 
      mu[i, t] <- beta0 + b1[i] * times[i, t]
      Mi[i, t] <- mu[i, t] + rho * (y[i, t - 1] - mu[i, t - 1])
    }
  }


  ## PRIORS
  # lambda

  lambda0 <- lambda * (1 - rho^2)

  lambda <- 1 / (sigma * sigma)
  sigma ~ dunif(0, 10)

  # rho

  rho ~ dunif(-1, 1)
  
  
  # beta0

  beta0 ~ dnorm(0, 0.01)

  # b1 prior

  for (i in 1: n){
    b1[i] ~ dnorm(0, lambda1)
  }
  lambda1 <- 1 / (sigma1 * sigma1)
  sigma1 ~ dunif(0, 10)
  

}"




#-------------------------------------------------------------------------------
# fit models
#-------------------------------------------------------------------------------

set.seed(1149989)
start_time <- Sys.time()

app_evidence_bridge_sampling <- c()

for(i in 1:5){

jags_H3 <- jags(data = list(y = y, n = n, sizes = sizes, times = times),
                parameters.to.save = c("beta0", "sigma", "rho", "sigma1", "b1"),
                model.file = textConnection(code_H3), n.chains = 1,
                n.iter = 1010000, n.burnin = 10000, n.thin = 100)

#-------------------------------------------------------------------------------
# specify unnormalized log posterior functions
#-------------------------------------------------------------------------------

library(mnormt)

log_posterior_H3 <- function(pars, data) {
  
  n <- length(data$y[, 1])
  T <- length(data$y[1, ])
  
  sizes <- c()
  for(i in 1 : n){
    sizes[i] <- T - sum(is.na(data$y[i,]))
  }
  
  beta0 <- pars["beta0"]           # extract parameter
  rho <- pars["rho"]               # extract parameter
  sigma <- pars["sigma"]           # extract parameter
  sigma1 <- pars["sigma1"]         # extract parameter
  b1 <- c()
  for(i in 1 : n){
    b1[i] <- pars[paste("b1[",i,"]", sep = '')]   # extract parameter
  }
  
  # likelihood
  log_likelihood <- 0
  for(i in 1 : n){
    X <- matrix(rep(1, sizes[i]), nrow = 1)
    Z <- matrix(c(data$times[i, 1 : sizes[i]]), nrow = 1)
    Sigma <- matrix(data = NA, nrow = sizes[i], ncol = sizes[i])
    for(t in 1 : sizes[i]){
      for(s in 1: sizes[i]){
        Sigma[t, s] <- rho^abs(t - s) * sigma^2 / (1 - rho^2)
      }
    }
    
    log_likelihood <- log_likelihood + dmnorm(x = data$y[i, 1 : sizes[i]], 
                                              mean = as.vector(t(X) %*% beta0) + as.vector(t(Z) %*% b1[i]), 
                                              varcov = Sigma, 
                                              log = TRUE) + # likelihood
      dnorm(b1[i], 0, sd = sigma1, log = TRUE) # random effects
  }
  
  
  out <- dnorm(beta0, mean = 0, sd = 10, log = TRUE) +         # prior
    dunif(sigma, 0, 10, log = TRUE) +                          # prior
    dunif(sigma1, 0, 10, log = TRUE) +                         # prior
    dunif(rho, 0, 10, log = TRUE) +                            # prior
    log_likelihood                                             # likelihood
  
  
  return(out)
  
}

#-------------------------------------------------------------------------------
# specify lower and upper bounds for the parameters
#-------------------------------------------------------------------------------

lb_H3 <- rep(-Inf, 14)
ub_H3 <- rep(Inf, 14)
names(lb_H3) <- names(ub_H3) <- c("beta0", "sigma", "sigma1", "rho",
                                  "b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", 
                                  "b1[6]", "b1[7]", "b1[8]", "b1[9]", "b1[10]")
lb_H3[["rho"]] <- -1
ub_H3[["rho"]] <- 1
lb_H3[["sigma"]] <- 0
ub_H3[["sigma"]] <- 10
lb_H3[["sigma1"]] <- 0
ub_H3[["sigma1"]] <- 10
#-------------------------------------------------------------------------------
# compute log marginal likelihoods
#-------------------------------------------------------------------------------

bridge_M3 <- bridge_sampler(samples = jags_H3,
                            log_posterior = log_posterior_H3,
                            data = list(y = y, n = n, sizes = sizes, times = times),
                            lb = lb_H3, ub = ub_H3, method = 'warp3',
                            repetitions = 1)
app_evidence_bridge_sampling <- c(app_evidence_bridge_sampling, bridge_M3$logml)
}
end_time <- Sys.time()

time_bs_M3 <- end_time - start_time

# summary(bridge_M3)

app_evidence_bridge_sampling_M3 <- app_evidence_bridge_sampling


save(app_evidence_bridge_sampling_M3, time_bs_M3, file = './bridgesample_example/results_bridgesample_M3.RData')
