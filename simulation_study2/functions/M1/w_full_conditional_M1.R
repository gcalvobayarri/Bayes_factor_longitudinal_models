w_full_conditional_M1 <- function(yy, tt, w_previous, bb1, beta00, ssigmau, 
                                  ssigma, rrho, ii, sizee, tau = 1){
  # sizee: number of individual measurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  w_posterior <- rep(NA, dim(yy)[2])
  
  w_posterior[1] <- rnorm(1, 
            mean = (ssigmau^2 * tau *
            (yy[ii, 1] - (beta00 + bb1[ii] * tt[ii, 1])) + 
            ssigma^2 * rrho * w_previous[2]) / (ssigma^2 + tau * ssigmau^2),
            sd = sqrt( (ssigma^2 * ssigmau^2) / (ssigma^2 + tau * ssigmau^2)  ))
  
  
  for(t in 2 : (sizee - 1)){
    w_posterior[t] <- rnorm(1, 
             mean = (ssigmau^2 * tau * 
             (yy[ii, t] - (beta00 + bb1[ii] * tt[ii, t])) +
             ssigma^2 * rrho * (w_posterior[t - 1] + w_previous[t + 1])) / 
             (ssigma^2 * (1 + rrho^2) + tau * ssigmau^2),
             sd = sqrt( (ssigma^2 * ssigmau^2) / (ssigma^2 * (1 + rrho^2) + tau * ssigmau^2)  ))
  }
  
  w_posterior[sizee] <- rnorm(1,
            mean = (ssigmau^2 * tau * 
            (yy[ii, sizee] - (beta00 + bb1[ii] * tt[ii, sizee])) + 
            ssigma^2 * rrho * w_posterior[sizee - 1]) /
            (ssigma^2 + tau * ssigmau^2),
            sd = sqrt( (ssigma^2 * ssigmau^2) / (ssigma^2 + tau * ssigmau^2)  ))
  
  return(w_posterior)
}
# wan <- rep(NA, dim(Y)[2])
# wan[1 : n[1]] <- c(rep(0, n[1]))
# wan <- w_full_conditional_M1(Y, times, wan, b1, beta0, sigmau, sigma, rho, 1, n[1]); print(wan)
# comprobar