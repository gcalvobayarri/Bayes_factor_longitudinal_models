sigma1_posterior_uniform_prior_random_trend_model <- function(sigma1_anterior, bb1, aa, bb) {
  
  nn <- length(bb1)
  
  # Siguiendo 2004_Robert_Casella
  sigma1_proposal <- runif(1, max(sigma1_anterior - 0.05, aa), min(sigma1_anterior + 0.05, bb))
  
  log_f_proposal <-  nn * (-log(sigma1_proposal)) + (- sum(bb1^2) / (2 * sigma1_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn * (-log(sigma1_anterior)) + (- sum(bb1^2) / (2 * sigma1_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma1_anterior, max(sigma1_proposal - 0.05, aa), min(sigma1_proposal + 0.05, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma1_proposal, max(sigma1_anterior - 0.05, aa), min(sigma1_anterior + 0.05, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigma1 <- sigma1_proposal
  }else{
    ssigma1 <- sigma1_anterior
  }
  
  return(ssigma1)
  
}

# load('./simulation study/data.RData')
# anterior <- 1
# anterior <- sigma1_posterior_uniform_prior_random_trend_model(anterior, c(-0.254, -0.1577, -0.04493, 0.36184, .21704), 0, 10)
# anterior

# comprobar sd de la posterior