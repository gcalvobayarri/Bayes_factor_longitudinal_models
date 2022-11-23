sigma_power_posterior_uniform_prior_random_trend_model <- function(yy, sigma_anterior, beta00, bb1, aa, bb, tau) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  diff_square <- matrix(data = NA, nrow = nn, ncol = TT)
  
  for(i in 1 : nn){
    for( t in 1 : TT){
      diff_square[i, t] <- (yy[i, t] - (beta00 + bb1[i] * (t - 1)))^2
    }
  }
  
  # Siguiendo 2004_Robert_Casella
  sigma_proposal <- runif(1, max(sigma_anterior - 0.05, aa), min(sigma_anterior + 0.05, bb))
  
  log_f_proposal <-  nn * TT * tau * (-log(sigma_proposal)) + (- tau * sum(diff_square) / (2 * sigma_proposal^2)  ) + log(1 / (bb - aa))
  
  log_f_anterior <-  nn * TT * tau * (-log(sigma_anterior)) + (- tau * sum(diff_square) / (2 * sigma_anterior^2)  ) + log(1 / (bb - aa))
  
  q_anterior_given_proposal <- dunif(sigma_anterior, max(sigma_proposal - 0.05, aa), min(sigma_proposal + 0.05, bb), log = T)
  
  q_proposal_given_anterior <- dunif(sigma_proposal, max(sigma_anterior - 0.05, aa), min(sigma_anterior + 0.05, bb), log = T)
  
  log_rho <- (log_f_proposal - log_f_anterior) + (q_anterior_given_proposal - q_proposal_given_anterior)
  
  rhoo <- min(1, exp(log_rho))
  
  if(runif(1) < rhoo){
    ssigma <- sigma_proposal
  }else{
    ssigma <- sigma_anterior
  }
  
  return(ssigma)
}

# load('./simulation study/data.RData')
# anterior <- 1
# anterior <- sigma_power_posterior_uniform_prior_random_trend_model(data, anterior, 1.97458, c(-0.254, -0.1577, -0.04493, 0.36184, .21704), 0, 10, 0)
# anterior #(tengo dudas)