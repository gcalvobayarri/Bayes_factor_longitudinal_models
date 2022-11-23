beta0_power_posterior_random_trend_model <- function(yy, bb1, ssigma, AA, tau) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  sum_b1t <- matrix(data = NA, nrow = nn, ncol = TT)
  for(i in 1 : nn){
    for(t in 0 : (TT - 1)){
      sum_b1t[i, t + 1] <- bb1[i] * t
    }
  }
  
  beta00 <- rnorm(1, AA^2 *  (nn * TT * tau * mean(yy) - tau * sum(sum_b1t)) / (tau * AA^2 * nn * TT + ssigma^2), 
                  sd = sqrt(  AA^2 * ssigma^2 / (tau * AA^2 * nn * TT + ssigma^2) ))
  
  return(beta00)
}

# load('./simulation study/data.RData')
# beta0_power_posterior_random_trend_model(data, c(-0.254, -0.1577, -0.04493, 0.36184, .21704), 1.10785, 10, 0)