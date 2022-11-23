beta0_power_posterior <- function(yy, bb0, bb1, ssigma, tau, AA) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  sum_b1_t <- 0
  for(i in 1 : length(bb1)){
    for(t in 0 : (TT -  1)){
      sum_b1_t <- sum_b1_t + bb1[i] * t
    }
  }
  
  beta00 <- rnorm(1, AA^2 * tau * (-TT * sum(bb0) - sum_b1_t + nn * TT * mean(yy)) / (AA^2 * tau * nn * TT + ssigma^2), 
                  sd = sqrt(  AA^2 * ssigma^2 / (10^2 * tau * nn * TT + ssigma^2) ))
  
  return(beta00)
}

# beta0_power_posterior(y, c(-1.25, 2.4, -10, 2, 2.2, 1.9, 4, -2.6, 1.5, 1, 2.6, -0.5),
# c(-0.025, 0.03, 0.09, 0.01, -0.0375, -0.005, -0.055, -0.0125, 0, -0.076, 0, 0.08), 0.541, 0)