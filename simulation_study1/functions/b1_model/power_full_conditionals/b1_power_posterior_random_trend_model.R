b1_power_posterior_random_trend_model <- function(yy, beta00, ssigma1, ssigma, tau) {
  
  bb1 <- c()
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  sum_t2 <- sum((1:(TT - 1))^2)
  
  
  
  for(i in 1 : dim(yy)[1]){
    bb1[i] <- rnorm(1, ssigma1^2 * tau * (-TT * (TT - 1) / 2 * beta00 + sum((0 : (TT - 1)) * yy[i, ])) / 
                      (ssigma1^2 * tau * sum_t2 + ssigma^2), 
                    sd = sqrt(  ssigma1^2 * ssigma^2 / (ssigma1^2 * tau * sum_t2 + ssigma^2) ))
  }
  
  
  
  return(bb1)
}

# load('./simulation study/data.RData')
# b1_power_posterior_random_trend_model(data, 1.97458, 0.34338, 1.10785, 0)