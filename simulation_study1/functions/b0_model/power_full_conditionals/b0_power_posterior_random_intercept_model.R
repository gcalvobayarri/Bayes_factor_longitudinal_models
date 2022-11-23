b0_power_posterior_random_intercept_model <- function(yy, beta00, ssigma0, ssigma, ttau) {
  
  bb0 <- c()
  TT <- dim(yy)[2]
  
  for(i in 1 : dim(yy)[1]){
    bb0[i] <- rnorm(1, ssigma0^2 * ttau * (-TT * beta00 + TT * mean(yy[i, ])) / (ssigma0^2 * ttau * TT + ssigma^2), 
                    sd = sqrt(  ssigma0^2 * ssigma^2 / (ssigma0^2 * ttau * TT + ssigma^2) ))
  }
  
  
  
  return(bb0)
}

# load('./simulation study/data.RData')
# b0_power_posterior_random_intercept_model(data, 2, 1.5, 0.5, 0)