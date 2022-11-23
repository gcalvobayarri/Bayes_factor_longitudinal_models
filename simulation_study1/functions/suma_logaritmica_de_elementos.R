suma_logaritmica_de_elementos <- function(a, n) {
  # a es vector de logaritmos
 
  log_sum <- a[1]
  for (i in 2 : n) {
    if(log_sum > a[i]){log_sum <- log_sum + log(1 + exp(a[i] - log_sum))}
    else{log_sum <- a[i] + log(1 + exp( log_sum - a[i]))}
  }
  
  return(log_sum)
}

# t <- proc.time()
# a[1] <- log(6.25)
# a[2] <- log(8.5)
# a[3] <- log(10)
# 
# suma_logaritmica_de_elementos(a, 3) == log(6.25 + 8.5 +10)
# tiempo_de_ejecucion <- proc.time() - t

# t <- proc.time()
# a <- log(rnorm(60000000, mean = 80, sd=1))
# 
# suma_logaritmica_de_elementos(a, 60000000)
# tiempo_de_ejecucion <- proc.time() - t


# a <- c(log(6), log(8))
# suma_logaritmica_de_elementos(a, 2) == log(6 + 8)