b1_power_posterior <- function(yy, beta00, bb0, ssigma1, ssigma, tau) {
  
  bb1 <- c()
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  sum_square_t <- sum((0 : (TT - 1))^2)
  
  sum_y_t <- rep(0, TT)
  for(i in 1 : nn){
    for(t in 0 : (TT -  1)){
      sum_y_t[i] <- sum_y_t[i] + yy[i, t + 1] * t
    }
  }
  
  for(i in 1 : nn){
    bb1[i] <- rnorm(1, ssigma1^2 * tau * (-TT * (TT - 1) / 2 * beta00 - bb0[i] * TT * (TT - 1) / 2 + sum_y_t[i]) / (ssigma1^2 * tau * sum_square_t + ssigma^2), 
                    sd = sqrt(  ssigma1^2 * ssigma^2 / (ssigma1^2 * tau * sum_square_t + ssigma^2) ))
  }
  
  
  
  return(bb1)
}

# b1_power_posterior(y, 8.098,  c(-1.25, 2.4, -10, 2, 2.2, 1.9, 4, -2.6, 1.5, 1, 2.6, -0.5), 0.054, 0.541, 0)