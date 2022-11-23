sigma0_posterior_uniform_prior <- function(sigma0_anterior, bb0, aa, bb) {
   
  nn <- length(bb0)
  
  # sigma0 
  sigma0_proposal <- runif(1, max(c(sigma0_anterior -.5,aa)), min(sigma0_anterior + .5, bb))
  rho <- 0
  for(n in 1 : nn){
    rho <- rho +     (log(1 / (sigma0_proposal * sqrt(2 * pi))) +
      (-0.5 * (bb0[n])^2 / sigma0_proposal^2)) -
    (log(1 / (sigma0_anterior * sqrt(2 * pi))) +
      (-0.5 * (bb0[n])^2 / sigma0_anterior^2))
      
      #(log(dnorm(bb0[n], 0, sd = sigma0_proposal))) - (log(dnorm(bb0[n], 0, sd = sigma0_anterior)))
  

    
    }
  rho <- rho + log(dunif(sigma0_anterior, max(c(sigma0_proposal -.5,aa)), min(sigma0_proposal +.5, bb))) - log(dunif(sigma0_proposal, max(c(sigma0_anterior -.5, aa)), min(sigma0_anterior +.5, bb)))
  
  ssigma0 <- sigma0_anterior + (sigma0_proposal - sigma0_anterior) * (runif(1)<exp(rho))
  
  return(ssigma0)
}

# anterior <- 1
# anterior <- sigma0_posterior_uniform_prior(anterior, c(-1.25, 2.4, -10, 2, 2.2, 1.9, 4, -2.6, 1.5, 1, 2.6, -0.5), 0, 10)
# anterior

# comprobar sd de la posterior
