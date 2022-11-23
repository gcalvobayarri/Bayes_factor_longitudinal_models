beta0_full_conditional_M2 <- function(yy, tt, bb1, ssigma, AA, 
                                      sizess, tau = 1){
  # sizess: number of individual measurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  err_matrix <- yy - ((t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt)
  
  
  beta00 <- rnorm(1, 
                  mean = AA^2 * tau * sum(err_matrix, na.rm = T) / (AA^2 * tau * sum(sizess) + ssigma^2),
                  sd = sqrt(AA^2 * ssigma^2 / (AA^2 * tau * sum(sizess) + ssigma^2) ))
  
  return(beta00)
}
# beta0_full_conditional_M2(Y, times, b1, sigma, 1, n)