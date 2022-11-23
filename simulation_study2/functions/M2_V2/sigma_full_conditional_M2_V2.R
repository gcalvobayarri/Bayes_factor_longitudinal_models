sigma_full_conditional_M2_V2 <- function(yy, tt, sigma_anterior, 
                                      beta00, bb0, bb1, aa, bb, sizess, tau = 1){
  # sizess: number of individual meassurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  sigma_proposal <- runif(1, max(sigma_anterior - 0.1, aa), min(sigma_anterior + 0.1, bb))
  
  sum_errors <- sum(
    (yy - (beta00 + (t(t(bb0)) %*% rep(1, dim(yy)[2])) + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt))^2, na.rm = T)
  
  log_f_proposal <-  tau * sum(sizess) * (-log(sigma_proposal)) + 
    (-tau * sum_errors) / (2 * sigma_proposal^2) +log(1 / (bb - aa))
  
  log_f_previous <-  tau * sum(sizess) * (-log(sigma_anterior)) + 
    (-tau * sum_errors) / (2 * sigma_anterior^2) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma_anterior, 
                                     max(sigma_proposal - 0.1, aa), min(sigma_proposal + 0.1, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma_proposal, 
                                     max(sigma_anterior - 0.1, aa), min(sigma_anterior + 0.1, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_previous) + 
    (q_anterior_given_proposal - q_proposal_given_anterior)
  
  ssigma <- sigma_anterior + (sigma_proposal - sigma_anterior) * (runif(1)<exp(log_rho))
  
  return(ssigma)
  
}