sigma_power_posterior <- function(yy, beta00, bb0, bb1, aa, bb, tau) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  sum_square <- 0
  for(i in 1 : nn){
    for(t in 0 :  (TT - 1)){
      sum_square <- sum_square + (yy[i, t + 1] - (beta00 + bb0[i] + bb1[i] * t))^2
    }
  }
  
  ssigma <- sqrt(rinvgamma(1,  shape = (nn * tau * TT) / 2 + aa, 
                           rate = (bb + 0.5 * tau * sum_square) ))
  
  return(ssigma)
}


# sigma_power_posterior(y, 8.098,  c(-1.25, 2.4, -10, 2, 2.2, 1.9, 4, -2.6, 1.5, 1, 2.6, -0.5),
# c(-0.025, 0.03, 0.09, 0.01, -0.0375, -0.005, -0.055, -0.0125, 0, -0.076, 0, 0.08), 0.01, 0.01, 0)