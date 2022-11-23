# 0. functions-------------------

source("./simulation_study2/functions/M2_V2/b0_full_conditional_M2_V2.R")
source("./simulation_study2/functions/M2_V2/b1_full_conditional_M2_V2.R")
source("./simulation_study2/functions/M2_V2/beta0_full_conditional_M2_V2.R")
source("./simulation_study2/functions/M2_V2/sigma_full_conditional_M2_V2.R")
source("./simulation_study2/functions/M2_V2/sigma0_full_conditional_M2_V2.R")
source("./simulation_study2/functions/M2_V2/sigma1_full_conditional_M2_V2.R")

source("./simulation_study2/functions/M2_V2/log_likelihood_M2_V2.R")

# 1. data-----------------

load('./simulation_study2/data_simulation_study2.RData')


N <- dim(Y)[1] # number of individuals
J <- dim(Y)[2] # maximum number of measurements

sizes <- c()
for(i in 1 : N){
  sizes[i] <- J - sum(is.na(Y[i,]))
}

# 2. Power posterior ------------------
# Iter
n.t <- 199 # number of different temperatures
c <- 5
#nburning <- 2000
nburning <- c(rep(30000, n.t + 1))
#niter_powpos <- 2000
# niter_powpos <- c(rep(200000, 5), rep(50000, 10), rep(30000, 5)) # probar con más
niter_powpos <- c(rep(20000, n.t + 1))

niter <- nburning + niter_powpos # burnin + niter_powpos necesita tener en total unas 400000 o más



tau_vector <- ((0 : n.t)/n.t)^c

beta0 <- c()
b1 <- matrix(data = NA, nrow = N, ncol =  sum(niter))
b0 <- matrix(data = NA, nrow = N, ncol =  sum(niter))
sigma <- c()
sigma0 <- c()
sigma1 <- c()

# beta0    N(0, A*A)
A <- 10

# sigmas ~ U(a, b)    (también sigma1)
a <- 0
b <- 10

# rho ~ U(-1, 1)

# 3. iterations----------------------
app_evidence_power_posterior <- c()

for (h in 1 : 10){
  
  log_likelihood <- c()
  # algorithm
  for(nt in 0 : n.t){
    cat(paste(nt/(n.t+1) * 100, '%'), " \r")
    if(nt == 0){
      
      beta0[1] <- 1
      #w[, , 1] <- matrix(data = rnorm(N * Y, 0, 3), nrow = N, ncol = Y)
      b0[, 1] <- runif(dim(Y)[1], 0, 10)
      b1[, 1] <- runif(dim(Y)[1], 0, 10)
      sigma[1] <- 1
      sigma1[1] <- 1
      sigma0[1] <- 1
      
      # likelihood
      log_likelihood[1] <- log_likelihood_M2_V2(Y, times, beta0[1], b0[, 1], 
                                                b1[, 1], sigma[1], sizes)
      
      
    }else{
      #revisar (hacer lo mismo que abajo)
      beta0[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(beta0[(sum(niter[1 : (nt)]) - niter[nt] + 1) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      sigma[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(sigma[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      sigma0[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(sigma0[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      sigma1[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(sigma1[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      
      for(j in 1 : N){
        b0[j, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
          mean(b0[j, (sum(niter[1 : (nt)]) - niter[nt] + 1) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
        
        b1[j, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
          mean(b1[j, (sum(niter[1 : (nt)]) - niter[nt] + 1) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
        
      }
      
      log_likelihood[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        log_likelihood_M2_V2(Y, times,
                          beta0[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1], 
                          b0[, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1],
                          b1[, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1],
                          sigma[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1], 
                          sizes)
      
      
    }
    
    
    for(i in 2 : niter[nt + 1]){
      ii <- sum(niter[1 : (nt + 1)]) - niter[nt + 1] + i # current iteration
      
      # beta0
      beta0[ii] <- 
        beta0_full_conditional_M2_V2(Y, times, bb0 = b0[, ii - 1],
                                  bb1 = b1[, ii - 1], 
                                  ssigma = sigma[ ii - 1], 
                                  AA = A, sizess = sizes, 
                                  tau = tau_vector[nt + 1])
      
      # sigma
      sigma[ii] <- 
        sigma_full_conditional_M2_V2(Y, times, 
                                  sigma[ii - 1], 
                                  beta0[ii],  b0[, ii - 1],
                                  b1[, ii - 1], a, b, sizes, 
                                  tau_vector[nt + 1])
      
      # sigma0
      sigma0[ii] <- 
        sigma0_full_conditional_M2_V2(sigma0[ii - 1], 
                                      b0[, ii - 1], a, b)
      
      # sigma1
      sigma1[ii] <- 
        sigma1_full_conditional_M2_V2(sigma1[ii - 1], 
                                   b1[, ii - 1], a, b)
      
      # b0   por aquí....
      b0[, ii] <-  
        b0_full_conditional_M2_V2(Y, times, 
                                  beta0[ii],
                                  b1[, ii - 1],
                                  sigma0[ii], 
                                  sigma[ii], sizes,
                                  tau_vector[nt + 1])
      
      # b1
      b1[, ii] <-  
        b1_full_conditional_M2_V2(Y, times, 
                               beta0[ii], b0[, ii],
                               sigma1[ii], 
                               sigma[ii],
                               tau_vector[nt + 1])
      
      
      log_likelihood[ii] <- 
        log_likelihood_M2_V2(yy=Y, tt = times,
                          beta00 = beta0[ii], bb0 = b0[, ii],
                          bb1 = b1[, ii],
                          ssigma = sigma[ii],
                          sizess = sizes)
    }
  } 
  # plot(beta0[400000:sum(niter)], type = 'l')
  # plot(sigma[100000:sum(niter)], type = 'l')
  # plot(sigmau[100000:sum(niter)], type = 'l')
  # plot(sigma1[100000:sum(niter)], type = 'l')
  
  expectation_loglikelihood_power_posterior <- c()
  
  for(i in 1 : length(tau_vector)){
    
    expectation_loglikelihood_power_posterior[i] <- 
      mean(log_likelihood[
        (sum(niter[1 : (i)]) - niter[i] + 1) : (sum(niter[1 : (i)]))]) 
    # hay algo mal...
  }
  
  
  # 6. Approximating evidence-----------------------
  log_evidence <- 0
  for(qq in 1 : (length(tau_vector) - 1)){
    log_evidence <- log_evidence + (tau_vector[qq + 1] - tau_vector[qq]) * 
      (expectation_loglikelihood_power_posterior[qq + 1] + 
         expectation_loglikelihood_power_posterior[qq]) / 2
  }
  
  
  
  app_evidence_power_posterior <- c(app_evidence_power_posterior, log_evidence)
  
  
  
  print(log_evidence)
  print(h)
}  

app_evidence_power_posterior_M2_V2 <- app_evidence_power_posterior

save(app_evidence_power_posterior_M2_V2, 
     file = "./simulation_study2/app_evidence_power_posterior_M2_V2.RData")
