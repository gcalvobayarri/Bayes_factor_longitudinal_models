beta0_power_posterior_random_intercept_model <- function(yy, bb0, ssigma, AA, ttau) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  beta00 <- rnorm(1, AA^2 * ttau * nn * TT * (-mean(bb0) + mean(yy)) / (AA^2 * ttau * nn * TT + ssigma^2), 
                  sd = sqrt(  AA^2 * ssigma^2 / (AA^2 * ttau * nn * TT + ssigma^2) ))
  
  return(beta00)
}

# load('./simulation study/data.RData')
# beta0_power_posterior_random_intercept_model(data, c(-1.58, -1.52, -0.53, 2.41, 1.44), 0.5, 10, 0)
#  