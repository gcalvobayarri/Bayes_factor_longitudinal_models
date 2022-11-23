sigma0_power_posterior_uniform_prior_random_intercept_model <- function(sigma0_anterior, bb0, aa, bb) {
  
  nn <- length(bb0)
  
  # Siguiendo 2004_Robert_Casella
  sigma0_proposal <- runif(1, max(sigma0_anterior - 0.05, aa), min(sigma0_anterior + 0.05, bb))
  
  log_f_proposal <-  nn * (-log(sigma0_proposal)) + (- sum(bb0^2) / (2 * sigma0_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn * (-log(sigma0_anterior)) + (- sum(bb0^2) / (2 * sigma0_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma0_anterior, max(sigma0_proposal - 0.05, aa), min(sigma0_proposal + 0.05, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma0_proposal, max(sigma0_anterior - 0.05, aa), min(sigma0_anterior + 0.05, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigma0 <- sigma0_proposal
  }else{
    ssigma0 <- sigma0_anterior
  }
  
  return(ssigma0)
  
}

# load('./simulation study/data.RData')
# anterior <- 1
# anterior <- sigma0_power_posterior_uniform_prior_random_intercept_model(anterior, c(-1.58, -1.52, -0.53, 2.41, 1.44), 0, 10)
# anterior

# comprobar sd de la posterior