sigma1_full_conditional_M3 <- function(sigma1_anterior, bb1, aa, bb) {
  
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)
  
  nn <- length(bb1)
  
  sigma1_proposal <- runif(1, max(sigma1_anterior - 0.1, aa), min(sigma1_anterior + 0.1, bb))
  
  log_f_proposal <-  nn * (-log(sigma1_proposal)) + (- sum(bb1^2) / (2 * sigma1_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn * (-log(sigma1_anterior)) + (- sum(bb1^2) / (2 * sigma1_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma1_anterior, max(sigma1_proposal - 0.1, aa), min(sigma1_proposal + 0.1, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma1_proposal, max(sigma1_anterior - 0.1, aa), min(sigma1_anterior + 0.1, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigma1 <- sigma1_proposal
  }else{
    ssigma1 <- sigma1_anterior
  }
  
  return(ssigma1)
  
}
# sig <- .8
# sig <- sigma1_full_conditional_M3(sig, b1, 0, 10); print(sig)