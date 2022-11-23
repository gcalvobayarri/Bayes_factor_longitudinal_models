sigma1_full_conditional_M2_V2 <- function(sigma1_anterior, bb1, aa, bb) {
  
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  nn <- length(bb1)
  
  sigma1_proposal <- runif(1, max(sigma1_anterior - 0.1, aa), min(sigma1_anterior + 0.1, bb))
  
  log_f_proposal <-  nn * (-log(sigma1_proposal)) + (- sum(bb1^2) / (2 * sigma1_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn * (-log(sigma1_anterior)) + (- sum(bb1^2) / (2 * sigma1_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma1_anterior, max(sigma1_proposal - 0.1, aa), min(sigma1_proposal + 0.1, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma1_proposal, max(sigma1_anterior - 0.1, aa), min(sigma1_anterior + 0.1, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  ssigma1 <- sigma1_anterior + (sigma1_proposal - sigma1_anterior) * (runif(1)<exp(log_rho))
  
  return(ssigma1)
  
}