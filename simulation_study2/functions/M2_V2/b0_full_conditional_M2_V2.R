b0_full_conditional_M2_V2 <- function(yy, tt, beta00, bb1, ssigma0, 
                                      ssigma, sizess, tau = 1){
  # sizess: number of individual measurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  bb0 <- c()
  
  ty_matrix <- (yy - (beta00 + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt))
  for(i in 1 : dim(yy)[1]){
    bb0[i] <- rnorm(1,
                    mean = (ssigma0^2 * tau * (sum(ty_matrix[i,], na.rm = T))) /
                      (ssigma^2 + tau * ssigma0^2 * sizess[i]),
                    sd = sqrt((ssigma * ssigma0)^2 / 
                                (ssigma^2 + tau * ssigma0^2 * sizess[i])))
  }
  
  return(bb0)
}
#b0_full_conditional_M2_V2(Y, times, beta0, b1, sigma0, sigma, n)
