set.seed(1)
N <- 10 #number of individuals
sigmau <- 1.5
rho <- .6
sigma1 <- .5
beta0 <- 2
sigma <- 2

n <- sample(10:70, N, replace = T) #sizes
times <- matrix(data = NA, nrow = N, ncol = max(n))
for(i in 1 : N){
  times[i, 1 : n[i]] <- runif(n[i], 0, 20)
  times[i, ] <- times[i, (order(times[i, ]))]
}



w <- matrix(data = NA, nrow = N, ncol = max(n))
for(i in 1 : N){
  w[i, 1] <- rnorm(1, 0, sd = sqrt(sigmau^2 / (1 - rho^2)))
  for(j in 2 : n[i]){
    w[i, j] <- rho * w[i, j - 1] + rnorm(1, 0, sd = sigmau)
  }
}

b1 <- rnorm(N, 0, sigma1)

Y <- matrix(data = NA, nrow = N, ncol = max(n))
for(i in 1 : N){
  for(j in 1 : n[i]){
    Y[i, j] <- rnorm(1, beta0 + b1[i] * times[i, j] + w[i, j], sigma)
  }
  
}

Time <- na.omit(c(t(times)))
yvec <- na.omit(c(t(Y)))
df <- data.frame(id = rep(1:N, times = n), Y = yvec, Time = Time)

library(ggplot2)
ggplot(df, aes(Time, Y, group = id)) + 
  geom_point() +
  geom_line() +
  #   theme(legend.title = element_blank())+
  coord_cartesian()+
  theme(legend.position = "none")+
  theme_minimal(base_size = 20, base_line_size = 1.0)+
  xlab('t')

save(Y, times, file = "./bridgesample_example/long_data_controversial.RData")
