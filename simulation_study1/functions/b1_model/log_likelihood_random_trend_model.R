log_likelihood_random_trend_model <- function(yy, bbeta0, bb1, ssigma, nn, TT){
  
  ln_pdf <- matrix(data = NA, nrow = nn, ncol = TT)
  log_likelihood_individual <- c()
  
  for(n in 1 : nn){
    
    for(t in 0 : (TT - 1)){
      ln_pdf[n, t + 1] <- log(1 / (ssigma * sqrt(2 * pi))) +
        (-0.5 * (yy[n, t + 1] - (bbeta0 + bb1[n] * (t)))^2 / ssigma^2)
      
      
      
    }
    log_likelihood_individual[n] <- (sum(ln_pdf[n, ]))
    
  }
  log_likelihood <- sum(log_likelihood_individual)
  return(log_likelihood)
}
