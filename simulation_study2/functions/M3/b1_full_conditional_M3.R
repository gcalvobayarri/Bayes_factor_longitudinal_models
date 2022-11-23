b1_full_conditional_M3 <- function(yy, tt, beta00, ssigma1, 
                                   ssigma, rrho, tau = 1){
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  bb1 <- c()
  
  t_diff_matrix <- tt[, 2 : dim(yy)[2]] - rrho * tt[, 1 : (dim(yy)[2] - 1)]
  y_diff_matrix <- rrho * yy[, 1 : (dim(yy)[2] - 1)] - yy[, 2 : dim(yy)[2]] +
    (1 - rrho) * beta00
  
  for(i in 1 : dim(yy)[1]){
    bb1[i] <- rnorm(1,
        mean = (ssigma1^2 * tau * (tt[i, 1] * (yy[i, 1] - beta00) - 
          sum(t_diff_matrix[i,] * y_diff_matrix[i,], na.rm = T))) /
          (ssigma^2 + tau * ssigma1^2 * ( (1 - rrho^2) * tt[i, 1]^2 + 
          sum(t_diff_matrix[i,]^2, na.rm = T))),
        sd = sqrt((ssigma * ssigma1)^2 / 
          (ssigma^2 + tau * ssigma1^2 * ( (1 - rrho^2) * tt[i, 1]^2 + 
          sum(t_diff_matrix[i,]^2, na.rm = T)))))
  }
  
  return(bb1)
}
# b1_full_conditional_M3(Y, times, beta0, sigma1, sigma, rho) # comprobar