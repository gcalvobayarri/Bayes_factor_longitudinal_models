sigma_power_posterior_uniform_prior <- function(yy, sigma_anterior, beta00, bb0, bb1, aa, bb, tau) {
  
  nn <- dim(yy)[1]
  TT <- dim(yy)[2]
  
  # sigma
  sigma_proposal <- runif(1, max(c(sigma_anterior - 0.1, aa)), min(sigma_anterior + 0.1, bb))
  ln_pdf_proposal <- matrix(data = NA, nrow = nn, ncol = TT)
  ln_pdf_anterior <- matrix(data = NA, nrow = nn, ncol = TT)
  # 
  rho <- 0
  for(n in 1 : nn){
    
    for(t in 0 : (TT - 1)){
      ln_pdf_proposal[n, t + 1] <- tau * (log(1 / (sigma_proposal * sqrt(2 * pi))) +
        (-0.5 * (yy[n, t + 1] - (beta00 + bb0[n] + bb1[n] * t))^2 / sigma_proposal^2))
      
      #   log(dnorm(yy[n, t + 1],
      # beta00 + bb0[n] + bb1[n] * t,  sigma_proposal ))
      
      
      
      # escribirlo como suma
      ln_pdf_anterior[n, t + 1] <-  tau * (log(1 / (sigma_anterior * sqrt(2 * pi))) +
        (-0.5 * (yy[n, t + 1] - (beta00 + bb0[n] + bb1[n] * t))^2 / sigma_anterior^2))
      
      #   log(dnorm(yy[n, t + 1],
      # beta00 + bb0[n] + bb1[n] * t,  sigma_anterior ))
      # 
      
      
    }
    rho <- rho + (sum(ln_pdf_proposal[n, ])) - sum(ln_pdf_anterior[n, ])
    
    #   rho <- rho + (log(dmnorm(yy[n,], beta00 + bb0[n] + bb1[n] * (0:(TT-1)), sigma_proposal^2 * diag(TT)))) - (log(dmnorm(yy[n,],beta00 + bb0[n] + bb1[n] * (0:(TT-1)), sigma_anterior^2 * diag(TT))))
  }
  rho <- rho + log(dunif(sigma_anterior, max(c(sigma_proposal - .1, aa)), max(sigma_proposal + .1, bb))) - log(dunif(sigma_proposal, max(c(sigma_anterior - .1, aa)), min(sigma_anterior + .1, bb)))
  
  ssigma <- sigma_anterior + (sigma_proposal - sigma_anterior) * (runif(1)<exp(rho))
  
  
  return(ssigma)
}


# anterior <- 1
# anterior <- sigma_power_posterior_uniform_prior(y, anterior, 8.098,  c(-1.25, 2.4, -10, 2, 2.2, 1.9, 4, -2.6, 1.5, 1, 2.6, -0.5),
# c(-0.025, 0.03, 0.09, 0.01, -0.0375, -0.005, -0.055, -0.0125, 0, -0.076, 0, 0.08), 0, 10, 0)
# anterior