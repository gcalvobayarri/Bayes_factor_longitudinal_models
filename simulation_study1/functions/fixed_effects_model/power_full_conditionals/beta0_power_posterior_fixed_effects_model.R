beta0_power_posterior_fixed_effects_model <- function(yy, ssigma, AA, tau) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  beta00 <- rnorm(1, tau * AA^2 * nn * TT * mean(yy) / (tau * AA^2 * nn * TT + ssigma^2), 
                  sd = sqrt(  AA^2 * ssigma^2 / (tau * AA^2 * nn * TT + ssigma^2) ))
  
  return(beta00)
}

# load('./simulation study/data.RData')
# beta0_power_posterior_fixed_effects_model(data, 1.7126, 10, 0)
#  