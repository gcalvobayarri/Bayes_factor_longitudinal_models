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

code_H1 <-
  "model{
  for (i in 1 : n){
    for (t in 1 : sizes[i]){
      y[i, t] ~ dnorm(Mi[i, t], lambda) 
      Mi[i, t] <- beta0 + b1[i] * times[i, t] + w[i, t]
    }
  }

  
  for (i in 1: n){
    w[i, 1] ~ dnorm(0, lambdau0)
    for (t in 2 : sizes[i]){
      w[i, t] ~ dnorm(rho * w[i, t-1], lambdau)
    }
  }

  ## PRIORS
  # lambda
  
  lambda <- 1 / (sigma * sigma)
  sigma ~ dunif(0, 10)

  # lambdau
  
  lambdau0 <- lambdau * (1 - rho^2)

  lambdau <- 1 / (sigmau * sigmau)
  sigmau ~ dunif(0, 10)

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

set.seed(12345)
start_time <- Sys.time()

app_evidence_bridge_sampling <- c()

for(i in 1:5){

jags_H1 <- jags(data = list(y = y, n = n, sizes = sizes, times = times),
                parameters.to.save = c("beta0", "sigma", "rho", "sigmau", 
                                       "sigma1", "b1"),
                model.file = textConnection(code_H1), n.chains = 1,
                n.iter = 1010000, n.burnin = 10000, n.thin = 100)



#-------------------------------------------------------------------------------
# specify unnormalized log posterior functions
#-------------------------------------------------------------------------------

library(mnormt)

log_posterior_H1 <- function(pars, data) {
  
  n <- length(data$y[, 1])
  T <- length(data$y[1, ])
  
  sizes <- c()
  for(i in 1 : n){
    sizes[i] <- T - sum(is.na(data$y[i,]))
  }
  
  beta0 <- pars["beta0"]           # extract parameter
  rho <- pars["rho"]               # extract parameter
  sigma <- pars["sigma"]           # extract parameter
  sigmau <- pars["sigmau"]         # extract parameter
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
    Sigmau <- matrix(data = NA, nrow = sizes[i], ncol = sizes[i])
    for(t in 1 : sizes[i]){
      for(s in 1: sizes[i]){
        Sigmau[t, s] <- rho^abs(t - s) * sigmau^2 / (1 - rho^2)
      }
    }
    
    log_likelihood <- log_likelihood + dmnorm(x = data$y[i, 1 : sizes[i]], 
                                              mean = as.vector(t(X) %*% beta0) + as.vector(t(Z) %*% b1[i]), 
                                              varcov = Sigmau + diag(sigma^2, sizes[i]), 
                                              log = TRUE) + # likelihood
      dnorm(b1[i], 0, sd = sigma1, log = TRUE) # random effects
  }
  
  
  out <- dnorm(beta0, mean = 0, sd = 10, log = TRUE) +         # prior
    dunif(sigma, 0, 10, log = TRUE) +                          # prior
    dunif(sigmau, 0, 10, log = TRUE) +                         # prior
    dunif(sigma1, 0, 10, log = TRUE) +                         # prior
    dunif(rho, 0, 10, log = TRUE) +                            # prior
    log_likelihood                                             # likelihood
  
  
  return(out)
  
}


#-------------------------------------------------------------------------------
# specify lower and upper bounds for the parameters
#-------------------------------------------------------------------------------

lb_H1 <- rep(-Inf, 15)
ub_H1 <- rep(Inf, 15)
names(lb_H1) <- names(ub_H1) <- c("beta0", "sigma", "sigmau", "sigma1", "rho",
                                  "b1[1]", "b1[2]", "b1[3]", "b1[4]", "b1[5]", 
                                  "b1[6]", "b1[7]", "b1[8]", "b1[9]", "b1[10]")
lb_H1[["rho"]] <- -1
ub_H1[["rho"]] <- 1
lb_H1[["sigma"]] <- 0
ub_H1[["sigma"]] <- 10
lb_H1[["sigmau"]] <- 0
ub_H1[["sigmau"]] <- 10
lb_H1[["sigma1"]] <- 0
ub_H1[["sigma1"]] <- 10



#-------------------------------------------------------------------------------
# compute log marginal likelihoods
#-------------------------------------------------------------------------------

bridge_M1 <- bridge_sampler(samples = jags_H1,
                            log_posterior = log_posterior_H1,
                            data = list(y = y, n = n, sizes = sizes, times = times),
                            lb = lb_H1, ub = ub_H1, method = 'warp3',
                            repetitions = 1)

app_evidence_bridge_sampling <- c(app_evidence_bridge_sampling, bridge_M1$logml)
}
end_time <- Sys.time()
time_bs_M1 <- end_time - start_time
# summary(bridge_M1)

app_evidence_bridge_sampling_M1 <- app_evidence_bridge_sampling



save(app_evidence_bridge_sampling_M1, time_bs_M1, file = './bridgesample_example/results_bridgesample_M1.RData')
