sigma_full_conditional_M1 <- function(yy, tt, sigma_anterior, 
                                      beta00, ww, bb1, aa, bb, sizess, tau = 1){
  # sizess: number of individual meassurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  sigma_proposal <- runif(1, max(sigma_anterior - 0.1, aa), min(sigma_anterior + 0.1, bb))
  
  sum_errors <- sum(
    (yy - (beta00 + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt + ww))^2, na.rm = T)
  
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
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigma <- sigma_proposal
  }else{
    ssigma <- sigma_anterior
  }
  
  return(ssigma)
  
}
# sig <- .2
# sig <- sigma_full_conditional_M1(Y, times, sig, beta0, w, b1, 0, 10, n); print(sig)