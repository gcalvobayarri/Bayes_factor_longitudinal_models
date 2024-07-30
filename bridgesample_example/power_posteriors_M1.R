# Power posterior

# 0. functions-------------------

source("./simulation_study2/functions/M1/w_full_conditional_M1.R")
source("./simulation_study2/functions/M1/b1_full_conditional_M1.R")
source("./simulation_study2/functions/M1/beta0_full_conditional_M1.R")
source("./simulation_study2/functions/M1/rho_full_conditional_M1.R")
source("./simulation_study2/functions/M1/sigma_full_conditional_M1.R")
source("./simulation_study2/functions/M1/sigma1_full_conditional_M1.R")
source("./simulation_study2/functions/M1/sigmau_full_conditional_M1.R")

source("./simulation_study2/functions/M1/log_likelihood_M1.R")

# 1. data-----------------

load("./bridgesample_example/long_data_controversial.RData")



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
#w <- array(data = NA, dim = c(N, Y, sum(niter))) # too large
w_actual <- matrix(data = NA, nrow = N, ncol = J)
w_anterior <- matrix(data = NA, nrow = N, ncol = J)
b1 <- matrix(data = NA, nrow = N, ncol =  sum(niter))
sigma <- c()
sigmau <- c()
sigma1 <- c()
rho <-  c()

# beta0    N(0, A*A)
A <- 10

# sigmas ~ U(a, b)    (también sigma1)
a <- 0
b <- 10

# rho ~ U(-1, 1)

# 3. iterations----------------------
app_evidence_power_posterior <- c()
set.seed(312665)
start_time_M1 <- Sys.time()
for (h in 1 : 5){
  
  log_likelihood <- c()
  # algorithm
  for(nt in 0 : n.t){
    cat(paste(nt/(n.t+1) * 100, '%'), " \r")
    if(nt == 0){
      
      beta0[1] <- 1
      #w[, , 1] <- matrix(data = rnorm(N * Y, 0, 3), nrow = N, ncol = Y)
      w_anterior <- matrix(data = rnorm(N * J, 0, 3), nrow = N, ncol = J)
      b1[, 1] <- runif(dim(Y)[1], 0, 10)
      sigma[1] <- 1
      sigmau[1] <- 1
      sigma1[1] <- 1
      rho[1] <- 0
      
      
      # likelihood
      log_likelihood[1] <- log_likelihood_M1(Y, times, beta0[1], b1[, 1], w_anterior, sigma[1], sizes)
      
      
    }else{
      #revisar (hacer lo mismo que abajo)
      beta0[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(beta0[(sum(niter[1 : (nt)]) - niter[nt] + 1) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      sigma[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(sigma[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      sigma1[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(sigma1[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      sigmau[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(sigmau[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      rho[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        mean(rho[(sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      # 
      # for(j in 1 : N){
      #   for(s in 1 : Y){
      #     w[j, s, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
      # mean(w[j, s, (sum(niter[1 : (nt)]) - niter[nt] + 1 ) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
      #   }
      # }
      for(j in 1 : N){
        b1[j, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
          mean(b1[j, (sum(niter[1 : (nt)]) - niter[nt] + 1) : (sum(niter[1 : (nt + 1)]) - niter[nt + 1])])
        
      }
      
      w_actual <- w_anterior
      
      log_likelihood[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1] <- 
        log_likelihood_M1(Y, times,
        beta0[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1], 
        b1[, sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1], 
        w_actual,
        sigma[sum(niter[1 : (nt + 1)]) - niter[nt + 1] + 1], 
        sizes)
      
      
    }

    
    for(i in 2 : niter[nt + 1]){
      ii <- sum(niter[1 : (nt + 1)]) - niter[nt + 1] + i # current iteration
      
      # beta0
      beta0[ii] <- 
        beta0_full_conditional_M1(Y, times, 
         bb1 = b1[, ii - 1], 
         ww = w_anterior, ssigma = sigma[ ii - 1], 
         AA = A, sizess = sizes, tau = tau_vector[nt + 1])
      # por aquí 16/12...
      # w
      for(n in 1 : N){
        w_actual[n, ] <- w_full_conditional_M1(Y, times, w_anterior[n,],
              b1[, ii - 1], 
              beta0[ii], 
              sigmau[ ii - 1], 
              sigma[ ii - 1], 
              rho[ ii - 1], n, 
              sizes[n], tau_vector[nt + 1])
        
      }
      
      # sigma
      sigma[ii] <- 
        sigma_full_conditional_M1(Y, times, 
          sigma[ii - 1], 
          beta0[ii], w_actual, 
          b1[, ii - 1], a, b, sizes, 
          tau_vector[nt + 1])
      
      # sigmau
      sigmau[ii] <- 
        sigmau_full_conditional_M1(sigmau[ii - 1], 
          w_actual, a, b, rho[ ii - 1], sizes)
      
      # sigma1
      sigma1[ii] <- 
        sigma1_full_conditional_M1(sigma1[ii - 1], 
         b1[, ii - 1], a, b)
      
      # b1
      b1[, ii] <-  
        b1_full_conditional_M1(Y, times, 
        beta0[ii], 
        sigma1[ii], 
        sigma[ii],
        w_actual, tau_vector[nt + 1])
      
      
      # rho
      rho[ii] <- 
        rho_full_conditional_M1(rho[ii - 1], 
         w_actual, sigmau[ii], sizes)
      
      w_anterior <- w_actual
      
      #likelihood el 500t mal...
      log_likelihood[ii] <- 
        log_likelihood_M1(yy=Y, tt = times,
         beta00 = beta0[ii],
         bb1 = b1[, ii], ww = w_actual,
         ssigma = sigma[ii],
         sizess = sizes)
    }
  } 
  # plot(beta0[500000:sum(niter)], type = 'l')
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
  
  
  # 4. Approximating evidence-----------------------
  log_evidence <- 0
  for(qq in 1 : (length(tau_vector) - 1)){
    log_evidence <- log_evidence + (tau_vector[qq + 1] - tau_vector[qq]) * 
      (expectation_loglikelihood_power_posterior[qq + 1] + 
      expectation_loglikelihood_power_posterior[qq]) / 2
  }
  
  # > log_evidence
  # [1] -172.8486
  
  app_evidence_power_posterior <- c(app_evidence_power_posterior, log_evidence)
  
  write.table(app_evidence_power_posterior, 
              file = "./bridgesample_example/app_evidence_power_posterior_M1.txt", 
              sep = ' ',
              row.names = F, col.names = F)

  print(log_evidence)
  print(h)
}  

app_evidence_power_posterior_M1 <- app_evidence_power_posterior
end_time_M1 <- Sys.time()
time_pp_M1 <- end_time_M1 - start_time_M1

save(app_evidence_power_posterior_M1, time_pp_M1,
     file = "./bridgesample_example/app_evidence_power_posterior_M1.RData")
