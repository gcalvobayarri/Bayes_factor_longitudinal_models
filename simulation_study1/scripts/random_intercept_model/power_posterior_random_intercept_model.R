# Power posterior random intercept model

# 0. functions-------------------

source("./simulation_study1/functions/b0_model/power_full_conditionals/b0_power_posterior_random_intercept_model.R")
source("./simulation_study1/functions/b0_model/power_full_conditionals/beta0_power_posterior_random_intercept_model.R")
source("./simulation_study1/functions/b0_model/power_full_conditionals/sigma_power_posterior_random_intercept_model.R")
source("./simulation_study1/functions/b0_model/power_full_conditionals/sigma0_power_posterior_random_intercept_model.R")

source("./simulation_study1/functions/b0_model/log_likelihood_random_intercept_model.R")



# 1. data-----------------

load('./simulation_study1/data_simulation_study1.RData')
y <- data


# 2. Power posterior ------------------
# prior


# Iter
library(LearnBayes)
library(invgamma)
library(DirichletReg)
nburning <- 20000
niter_powpos <- 30000
niter <- nburning + niter_powpos # burnin + niter_powpos

n.t <- 199 # number of different temperatures
c <- 5
tau_vector <- ((0 : n.t)/n.t)^c

beta0 <- c()
b0 <- matrix(data = NA, nrow = dim(y)[1], ncol =  (n.t + 1) * niter)
sigma <- c()
sigma0 <- c()

A <- 10

# sigmas ~ U(a, b)    (también sigma0)
a <- 0
b <- 10



# 3. iterations----------------------
app_evidence_power_posterior <- c()

for (h in 1 : 10) {
  
  # algorithm
  for(nt in 0 : n.t){
    print(nt)
    
    if(nt == 0){
      
      beta0[1] <- 1
      b0[, 1] <- runif(dim(y)[1], 0, 10)
      sigma[1] <- 1
      sigma0[1] <- 1

    }else{
      
      beta0[nt * niter + 1] <- mean(beta0[((nt - 1) * niter + 1 ) : (nt * niter)])
      sigma[nt * niter + 1] <- mean(sigma[((nt - 1) * niter + 1 ) : (nt * niter)])
      sigma0[nt * niter + 1] <- mean(sigma0[((nt - 1) * niter + 1 ) : (nt * niter)])

      for(j in 1 : dim(y)[1]){
        b0[j, nt * niter + 1] <- mean(b0[j, ((nt - 1) * niter + 1 ) : (nt * niter)])
      }
      
      
    }
    #niter <- 60000
    # nt <- 19
    
    for(i in 2 : niter){
      
      # beta0
      beta0[nt * niter + i] <- beta0_power_posterior_random_intercept_model(y, b0[, nt * niter + i - 1], sigma[ nt * niter + i - 1], tau_vector[nt + 1], A)
      
      # b0
      b0[, nt * niter + i] <-  b0_power_posterior_random_intercept_model(y, beta0[nt * niter + i], sigma0[nt * niter + i -1], sigma[nt * niter + i -1], tau_vector[nt + 1])
      
      # sigma
      sigma[nt * niter + i] <- sigma_power_posterior_uniform_prior_random_intercept_model(y, sigma[nt * niter + i - 1], beta0[nt * niter + i], b0[, nt * niter + i], a, b, tau_vector[nt + 1])
      
      # sigma0
      sigma0[nt * niter + i] <- sigma0_power_posterior_uniform_prior_random_intercept_model(sigma0[nt * niter + i - 1], b0[, nt * niter + i], a, b)
      
    }
    
  } 
  
  # plot(beta0, type = 'l')
  # plot(sigma, type = 'l')
  # plot(sigma0, type = 'l')
  # plot(sigma1, type = 'l')
  # 
  # plot(beta0[1:niter], type = 'l')
  # plot(sigma[1:niter], type = 'l')
  # plot(sigma0[1:niter], type = 'l')
  # plot(sigma1[1:niter], type = 'l')
  # 
  # plot(beta0[((n.t) * niter + 1):((n.t + 1) * niter)], type = 'l')
  # plot(sigma[((n.t) * niter + 1):((n.t + 1) * niter)], type = 'l')
  # plot(sigma0[((n.t) * niter + 1):((n.t + 1) * niter)], type = 'l')
  # plot(sigma1[((n.t) * niter + 1):((n.t + 1) * niter)], type = 'l') # hay algo mal
  
  # > mean(sigma1[119000:120000])
  # [1] 0.1464444
  # > sd(sigma1[119000:120000])
  # [1] 0.03254891
  # > mean(sigma0[119000:120000])
  # [1] 4.079715
  # > sd(sigma0[119000:120000])
  # [1] 0.9665205
  # > sd(sigma1[114001:120000])
  # [1] 0.03234646
  # > mean(sigma1[114001:120000])
  # [1] 0.1457734
  # > mean(beta0[119000:120000])
  # [1] 6.848038
  # > sigma(beta0[119000:120000])
  # Error: $ operator is invalid for atomic vectors
  # > sd(beta0[119000:120000])
  # [1] 0.2520037
  
  # 4. Likelihood (conditional)---------------------
  
  #ln_pdf <- array(NA, dim = c(dim(y)[1], dim(y)[2], (n.t + 1) * niter_powpos))
  log_likelihood_individual <- matrix(NA, nrow = dim(y)[1], ncol = (n.t + 1) * niter_powpos)
  log_likelihood <- c()
  for(nt in 0 : n.t){
    for(i in 1 : niter_powpos){
      # for(n in 1 : dim(y)[1]){
      #   
      #   for(t in 0 : (dim(y)[2] - 1)){
      #     ln_pdf[n, t + 1, nt * niter_powpos + i] <- log(1 / (sigma[nt * niter + nburning + i] * sqrt(2 * pi))) +
      #       (-0.5 * (y[n, t + 1] - (beta0[nt * niter + nburning + i] + b0[n, nt * niter + nburning + i] + b1[n, nt * niter + nburning + i] * t))^2 / sigma[nt * niter + nburning + i]^2)
      #     
      #     
      #     
      #   }
      #   log_likelihood_individual[n, nt * niter_powpos + i] <- (sum(ln_pdf[n, , nt * niter_powpos + i]))
      #   
      # }
      # log_likelihood[ nt * niter_powpos + i] <- sum(log_likelihood_individual[, nt * niter_powpos + i])
      
      
      log_likelihood[ nt * niter_powpos + i] <- log_likelihood_random_intercept_model(y, beta0[nt * niter + nburning + i], 
            b0[, nt * niter + nburning + i], sigma[nt * niter + nburning + i], dim(y)[1], dim(y)[2])
      
    }
  }
  
  # 5. Expectation likelihood from the power posteriors----------
  
  expectation_loglikelihood_power_posterior <- c()
  
  for(i in 1 : length(tau_vector)){
    
    expectation_loglikelihood_power_posterior[i] <- mean(log_likelihood[((i-1) * niter_powpos + 1) : (i * niter_powpos)])
    
  }
  
  
  # 6. Approximating evidence-----------------------
  log_evidence <- 0
  for(qq in 1 : (length(tau_vector) - 1)){
    log_evidence <- log_evidence + (tau_vector[qq + 1] - tau_vector[qq]) * (expectation_loglikelihood_power_posterior[qq + 1] + expectation_loglikelihood_power_posterior[qq]) / 2
  }
  
  print(log_evidence)
  app_evidence_power_posterior <- c(app_evidence_power_posterior, log_evidence)
  
}




app_evidence_power_posterior_random_intercept_model_simulated_data <- app_evidence_power_posterior

save(app_evidence_power_posterior_random_intercept_model_simulated_data, 
     file = "./simulation_study1/results/app_evidence_power_posterior_random_intercept_model_simulated_data.RData")
load("./simulation_study1/results/app_evidence_power_posterior_random_intercept_model_simulated_data.RData")


