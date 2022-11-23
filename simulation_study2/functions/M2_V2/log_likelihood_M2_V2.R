log_likelihood_M2_V2 <- function(yy, tt, beta00, bb0, bb1, ssigma, sizess){
  # sizess: number of individual measurements
  
  diff_t <- yy - (beta00 + (t(t(bb0)) %*% rep(1, dim(yy)[2])) + 
                    (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt)
  
  log_like <- sum(sizess) * (-log(ssigma * sqrt(2 * pi))) + 
    (- sum(diff_t^2, na.rm = T) / (2 * ssigma^2)) 
  
  return(log_like)
  
}