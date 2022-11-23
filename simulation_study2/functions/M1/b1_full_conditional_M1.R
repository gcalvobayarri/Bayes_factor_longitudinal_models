b1_full_conditional_M1 <- function(yy, tt, beta00, ssigma1, 
                                   ssigma, ww, tau = 1){
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  bb1 <- c()
  
  ty_matrix <- tt * (yy - (beta00 + ww))
  t2_matrix <- tt * tt
  
  for(i in 1 : dim(yy)[1]){
    bb1[i] <- rnorm(1,
          mean = (ssigma1^2 * tau * (sum(ty_matrix[i,], na.rm = T))) /
            (ssigma^2 + tau * ssigma1^2 * sum(t2_matrix[i,], na.rm = T)),
          sd = sqrt((ssigma * ssigma1)^2 / 
            (ssigma^2 + tau * ssigma1^2 * sum(t2_matrix[i,], na.rm = T))))
  }
  
  return(bb1)
}
# b1_full_conditional_M1(Y, times, beta0, sigma1, sigma, w) # comprobar