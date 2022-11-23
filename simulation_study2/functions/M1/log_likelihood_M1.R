log_likelihood_M1 <- function(yy, tt, beta00, bb1, ww, ssigma, sizess){
  # sizess: number of individual measurements
  
  diff_t <- yy - (beta00 + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt + ww)
  
  log_like <- sum(sizess) * (-log(ssigma * sqrt(2 * pi))) + 
    (- sum(diff_t^2, na.rm = T) / (2 * ssigma^2)) 
  
  return(log_like)
  
}
