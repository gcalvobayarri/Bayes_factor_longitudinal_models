rho_full_conditional_M3 <- function(yy, tt, rho_anterior, 
                                    beta00, ssigma, bb1, sizess, tau = 1){
  
  # sizess: number of individual meassurements
  # tau : temperature of the power posterior between 0 (prior) and 1 (posterior)

  rho_proposal <- runif(1, max(c(rho_anterior - 0.05, -1)), min(rho_anterior + 0.05, 1))
  
  error1_matrix <- 
    (yy - (beta00 + (t(t(bb1)) %*% rep(1, dim(yy)[2])) * tt))
  
  sum_error1 <- sum(error1_matrix[, 1]^2, na.rm = T)
  
  error_matrix_proposal <- 
    yy[, 2 : dim(yy)[2]] - (    (1 - rho_proposal) * beta00 + 
    (t(t(bb1)) %*% rep(1, dim(yy)[2] - 1)) * 
    (tt[, 2 : dim(yy)[2]] - rho_proposal * tt[, 1 : (dim(yy)[2] - 1)]) + 
    rho_proposal * yy[, 1 : (dim(yy)[2] - 1)]     )
  
  sum_error_proposal <- sum(error_matrix_proposal^2, na.rm = T)
  
  error_matrix_anterior <- 
    yy[, 2 : dim(yy)[2]] - (    (1 - rho_anterior) * beta00 + 
                                  (t(t(bb1)) %*% rep(1, dim(yy)[2] - 1)) * 
                                  (tt[, 2 : dim(yy)[2]] - rho_anterior * tt[, 1 : (dim(yy)[2] - 1)]) + 
                                  rho_anterior * yy[, 1 : (dim(yy)[2] - 1)]     )
  
  sum_error_anterior <- sum(error_matrix_anterior^2, na.rm = T)
  
  # Metropolis H (por aqui)
  
  log_f_proposal <- tau * (dim(yy)[1] / 2) * (-log(ssigma^2 / (1 - rho_proposal^2))) + 
    (-tau * (1 - rho_proposal^2) * sum_error1) / (2 * ssigma^2) +
    tau * sum(sizess - 1) * (-log(ssigma)) + 
    (-tau * sum_error_proposal) / (2 * ssigma^2) + log(1 / 2)
  
  log_f_anterior <-  tau * (dim(yy)[1] / 2) * (-log(ssigma^2 / (1 - rho_anterior^2))) + 
    (-tau * (1 - rho_anterior^2) * sum_error1) / (2 * ssigma^2) +
    tau * sum(sizess - 1) * (-log(ssigma)) + 
    (-tau * sum_error_anterior) / (2 * ssigma^2) + log(1 / 2)
  
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
# sig <- rho_full_conditional_M3(Y, times, sig, beta0, .5, b1, n); print(sig)