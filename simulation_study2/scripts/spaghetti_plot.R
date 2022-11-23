load('./simulation_study2/data_simulation_study2.RData')

N <- dim(Y)[1]
J <- dim(Y)[2]

sizes <- c()
for(i in 1 : N){
  sizes[i] <- J - sum(is.na(Y[i,]))
}

df <- data.frame(y = na.omit(c(t(Y))), time = na.omit(c(t(times))), id = rep(1 : N, times = sizes))

library(ggplot2)
library(tidyverse); library(scales)
ggplot(df, aes(time, y, group = id)) +
  geom_point() +
  geom_line() +
  coord_cartesian()+
  theme_minimal(base_size = 20, base_line_size = 1.0)+
  xlab('t')#+
  # scale_x_continuous(
  #   breaks = seq(from=0,to=19,by=4),
  #   minor_breaks = seq(from=1,to=20,by=4)
  #)
