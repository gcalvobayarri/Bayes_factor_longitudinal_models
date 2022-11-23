sigmau_full_conditional_M1 <- function(sigmau_anterior, ww, aa, bb, rrho, sizess) {
  
  # sizess: number of individual meassurements
  
  nn <- length(sizess)
  TT <- max(sizess)
  
  # some calculations
  rho_w <- matrix(data = NA, nrow = nn, ncol = TT - 1)
  for(i in 1 : nn){
    rho_w[i, ] <- rrho * ww[i, 1 : (TT - 1)]
  }
  
  
  # Metropolis H
  sigmau_proposal <- runif(1, max(sigmau_anterior - 0.1, aa), min(sigmau_anterior + 0.1, bb))
  
  log_f_proposal <- nn / 2 * (-log(sigmau_proposal^2 / (1 - rrho^2))) + 
    (- sum(ww[, 1]^2) / (2 * (sigmau_proposal^2 / (1 - rrho^2)))) +    
     (sum(sizess) - nn) * (-log(sigmau_proposal)) + 
    (- sum((ww[1 : nn, 2 : TT] - rho_w)^2, na.rm = T) / (2 * sigmau_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn / 2 * (-log(sigmau_anterior^2 / (1 - rrho^2))) + 
    (- sum(ww[, 1]^2) / (2 * (sigmau_anterior^2 / (1 - rrho^2)))) +    
    (sum(sizess) - nn) * (-log(sigmau_anterior)) + 
    (- sum((ww[1 : nn, 2 : TT] - rho_w)^2, na.rm = T) / (2 * sigmau_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigmau_anterior, max(sigmau_proposal - 0.1, aa), min(sigmau_proposal + 0.1, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigmau_proposal, max(sigmau_anterior - 0.1, aa), min(sigmau_anterior + 0.1, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigmau <- sigmau_proposal
  }else{
    ssigmau <- sigmau_anterior
  }
  
  return(ssigmau)
}
# sig <- .2
# sig <- sigmau_full_conditional_M1(sig, w, 0, 10, rho, n); print(sig)