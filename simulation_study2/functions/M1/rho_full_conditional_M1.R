rho_full_conditional_M1 <- function(rho_anterior, ww, ssigmau, sizess){
  
  # sizess: number of individual meassurements
  
  nn <- length(sizess)
  TT <- max(sizess)
  
  rho_proposal <- runif(1, max(c(rho_anterior - 0.05, -1)), min(rho_anterior + 0.05, 1))
  
  # some calculations (a ver si esto lo que está mal)
  rho_w_proposal <- matrix(data = NA, nrow = nn, ncol = TT - 1)
  for(i in 1 : nn){
    rho_w_proposal[i, ] <- rho_proposal * ww[i, 1 : (TT - 1)]
  }
  
  rho_w_anterior <- matrix(data = NA, nrow = nn, ncol = TT - 1)
  for(i in 1 : nn){
    rho_w_anterior[i, ] <- rho_anterior * ww[i, 1 : (TT - 1)]
  }
  
  # Metropolis H
  
  log_f_proposal <- nn / 2 * log(1 - rho_proposal^2)  + (- ((1 - rho_proposal^2) * sum( ww[, 1]^2) + sum((ww[1:nn,2:TT] - rho_w_proposal)^2, na.rm = T)) / (2 * ssigmau^2)  ) + log(1 / 2)
  
  log_f_anterior <-  nn / 2 * log(1 - rho_anterior^2)  + (- ((1 - rho_anterior^2) * sum( ww[, 1]^2) + sum((ww[1:nn,2:TT] - rho_w_anterior)^2, na.rm = T)) / (2 * ssigmau^2)  ) + log(1 / 2)
  
  q_anterior_given_proposal <- dunif(rho_anterior, max(rho_proposal - 0.05, -1), min(rho_proposal + 0.05, 1), log = T)
  
  q_proposal_given_anterior <- dunif(rho_proposal, max(rho_anterior - 0.05, -1), min(rho_anterior + 0.05, 1), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    rrho <- rho_proposal
  }else{
    rrho <- rho_anterior
  }
  
  return(rrho)
}
# sig <- .2
# sig <- rho_full_conditional_M1(sig, w, sigmau, n); print(sig)