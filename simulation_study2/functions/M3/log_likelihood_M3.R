log_likelihood_M3 <- function(yy, tt, beta00, bb1, rrho, ssigma, sizess){
  # sizess: number of individual measurements
  
  error1_matrix <- 
    (yy - (beta00 + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt))
  
  sum_error1 <- sum(error1_matrix[, 1]^2, na.rm = T)
  
  error_matrix <- 
    yy[, 2 : dim(yy)[2]] - (    (1 - rrho) * beta00 + 
                                  (t(t(bb1)) %*% rep(1, dim(yy)[2] - 1)) * 
                                  (tt[, 2 : dim(yy)[2]] - rrho * tt[, 1 : (dim(yy)[2] - 1)]) + 
                                  rrho * yy[, 1 : (dim(yy)[2] - 1)]     )
  
  sum_error <- sum(error_matrix^2, na.rm = T)
  
  log_like <- (dim(yy)[1] / 2) * (-log(ssigma^2 / (1 - rrho^2) * 2 * pi)) + 
    (-(1 - rrho^2) * sum_error1) / (2 * ssigma^2) +
     sum(sizess - 1) * (-log(ssigma * sqrt(2 * pi))) + 
    (-sum_error) / (2 * ssigma^2)
  
  return(log_like)
  
}