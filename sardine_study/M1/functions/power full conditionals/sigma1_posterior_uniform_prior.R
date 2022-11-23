sigma1_posterior_uniform_prior <- function(sigma1_anterior, bb1, aa, bb) {
  
  nn <- length(bb1)
  
  # sigma0 
  sigma1_proposal <- runif(1, max(c(sigma1_anterior - 0.05,aa)), min(sigma1_anterior + 0.05, bb))
  rho <- 0
  for(n in 1 : nn){
    rho <- rho + (log(1 / (sigma1_proposal * sqrt(2 * pi))) +
        (-0.5 * (bb1[n])^2 / sigma1_proposal^2)) -
      (log(1 / (sigma1_anterior * sqrt(2 * pi))) +
         (-0.5 * (bb1[n])^2 / sigma1_anterior^2))
      
      
      #(log(dnorm(bb1[n], 0, sd = sigma1_proposal))) - (log(dnorm(bb1[n], 0, sd = sigma1_anterior)))
  
    
    
    }
  rho <- rho + log(dunif(sigma1_anterior, max(c(sigma1_proposal - 0.05, aa)), min(sigma1_proposal + 0.05, bb))) - log(dunif(sigma1_proposal, max(c(sigma1_anterior - 0.05, aa)), min(sigma1_anterior + 0.05, bb)))
  
  ssigma1 <- sigma1_anterior + (sigma1_proposal - sigma1_anterior) * (runif(1)<exp(rho))
  
  return(ssigma1)
}

# anterior <- 1
# anterior <- sigma1_posterior_uniform_prior(anterior, c(-0.025, 0.03, 0.09, 0.01, -0.0375, -0.005, -0.055, -0.0125, 0, -0.076, 0, 0.08), 0, 10)
# anterior

# comprobar sd de la posterior