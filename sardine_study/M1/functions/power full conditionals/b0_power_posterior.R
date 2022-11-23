b0_power_posterior <- function(yy, beta00, bb1, ssigma0, ssigma, tau) {
  
  bb0 <- c()
  TT <- dim(yy)[2]
  
  sum_b1_t <- rep(0, TT)
  for(i in 1 : length(bb1)){
    for(t in 0 : (TT -  1)){
      sum_b1_t[i] <- sum_b1_t[i] + bb1[i] * t
    }
  }
  
  for(i in 1 : length(bb1)){
    bb0[i] <- rnorm(1, ssigma0^2 * tau * (-TT * beta00 - sum_b1_t[i] + TT * mean(yy[i, ])) / (ssigma0^2 * tau * TT + ssigma^2), 
                    sd = sqrt(  ssigma0^2 * ssigma^2 / (ssigma0^2 * tau * TT + ssigma^2) ))
  }
  
  
  
  return(bb0)
}

# b0_power_posterior(y, 8.098,  c(-0.025, 0.03, 0.09, 0.01, -0.0375, -0.005, -0.055, -0.0125, 0, -0.076, 0, 0.08), 4.234, 0.541, 0)
