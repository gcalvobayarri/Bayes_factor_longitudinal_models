beta0_full_conditional_M3 <- function(yy, tt, bb1, ssigma, rrho, AA, 
                                      sizess, tau = 1){
  # sizess: number of individual measurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  err_matrix <- yy - ((t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt)
  
  err_matrix2 <- rrho * yy[, 1 : (dim(yy)[2] - 1)] - yy[, 2 : dim(yy)[2]] + 
    ((t(t(bb1)) %*% rep(1, dim(yy)[2] - 1)) * 
    (tt[, 2 : dim(yy)[2]] - rrho * tt[, 1 : (dim(yy)[2] - 1)]))
  
  
  beta00 <- rnorm(1, 
    mean = AA^2 * tau * ((1 - rrho^2) *
    sum(err_matrix[, 1], na.rm = T) - (1 - rrho) * sum(err_matrix2, na.rm = T)) /
    (AA^2 * tau * ( (1 - rrho^2) * dim(yy)[1] + (1 - rrho)^2 * sum(sizess - 1)) + ssigma^2),
    sd = sqrt(AA^2 * ssigma^2 /
    (AA^2 * tau * ( (1 - rrho^2) * dim(yy)[1] + (1 - rrho)^2 * sum(sizess - 1)) + ssigma^2) ))
  
  return(beta00)
}
# beta0_full_conditional_M3(Y, times, b1, sigma, rho, 1, n)