load('./simulation_study1/data_simulation_study1.RData')
data[2,1]


df <- data.frame(y = as.numeric(as.vector(t(data))), time = rep(0:9, times = 5), id = rep(1 : 5, each = 10))

library(ggplot2)
library(tidyverse); library(scales)
ggplot(df, aes(time, y, group = id)) +
  geom_point() +
  geom_line() +
  coord_cartesian()+
  theme_minimal(base_size = 20, base_line_size = 1.0)+
  xlab('t')+
  scale_x_continuous(
    breaks = seq(from=0,to=9,by=2),
    minor_breaks = seq(from=1,to=10,by=2)
  )
