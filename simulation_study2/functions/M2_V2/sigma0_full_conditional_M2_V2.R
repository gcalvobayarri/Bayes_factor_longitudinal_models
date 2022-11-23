sigma0_full_conditional_M2_V2 <- function(sigma0_anterior, bb0, aa, bb) {
  
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  nn <- length(bb0)
  
  sigma0_proposal <- runif(1, max(sigma0_anterior - 0.1, aa), min(sigma0_anterior + 0.1, bb))
  
  log_f_proposal <-  nn * (-log(sigma0_proposal)) + (- sum(bb0^2) / (2 * sigma0_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn * (-log(sigma0_anterior)) + (- sum(bb0^2) / (2 * sigma0_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma0_anterior, max(sigma0_proposal - 0.1, aa), min(sigma0_proposal + 0.1, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma0_proposal, max(sigma0_anterior - 0.1, aa), min(sigma0_anterior + 0.1, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  ssigma0 <- sigma0_anterior + (sigma0_proposal - sigma0_anterior) * (runif(1)<exp(log_rho))
  
  return(ssigma0)
  
}
# sig <- 1
# for(i in 1 : 10000){sig <- sigma0_full_conditional_M2_V2(sig, b0, 0, 10)}
