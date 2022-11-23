sigma_full_conditional_M3 <- function(yy, tt, sigma_anterior, 
                                beta00, bb1, rrho, aa, bb, sizess, tau = 1){
  # sizess: number of individual meassurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  sigma_proposal <- 
    runif(1, max(sigma_anterior - 0.1, aa), min(sigma_anterior + 0.1, bb))
  
  error1_matrix <- 
    (yy - (beta00 + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt))
  
  sum_error1 <- sum(error1_matrix[, 1]^2, na.rm = T)
  
  error_matrix <- 
    yy[, 2 : dim(yy)[2]] - (    (1 - rrho) * beta00 + 
      (t(t(bb1)) %*% rep(1, dim(yy)[2] - 1)) * 
      (tt[, 2 : dim(yy)[2]] - rrho * tt[, 1 : (dim(yy)[2] - 1)]) + 
      rrho * yy[, 1 : (dim(yy)[2] - 1)]     )
  
  sum_error <- sum(error_matrix^2, na.rm = T)
  
  
  log_f_proposal <-  tau * (dim(yy)[1] / 2) * (-log(sigma_proposal^2 / (1 - rrho^2))) + 
    (-tau * (1 - rrho^2) * sum_error1) / (2 * sigma_proposal^2) +
    tau * sum(sizess - 1) * (-log(sigma_proposal)) + 
    (-tau * sum_error) / (2 * sigma_proposal^2) +log(1 / (bb - aa))
  
  log_f_previous <-  tau * (dim(yy)[1] / 2) * (-log(sigma_anterior^2 / (1 - rrho^2))) + 
    (-tau * (1 - rrho^2) * sum_error1) / (2 * sigma_anterior^2) +
    tau * sum(sizess - 1) * (-log(sigma_anterior)) + 
    (-tau * sum_error) / (2 * sigma_anterior^2) +log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma_anterior, 
                                     max(sigma_proposal - 0.1, aa), min(sigma_proposal + 0.1, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma_proposal, 
                                     max(sigma_anterior - 0.1, aa), min(sigma_anterior + 0.1, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_previous) + 
    (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigma <- sigma_proposal
  }else{
    ssigma <- sigma_anterior
  }
  
  return(ssigma)
  
}
# sig <- .2
# sig <- sigma_full_conditional_M3(Y, times, sig, beta0, b1, rho, 0, 10, n); print(sig)