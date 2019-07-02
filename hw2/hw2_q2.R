library(dplyr)

data <- read.csv("corruption.csv")
cutoff <- 10188
data_Q2 <- filter(data, pop >= cutoff - 1000, pop <= cutoff + 1000)
rdplot(y = data_Q2$fpm, x = data_Q2$pop,
       x.label = "population", y.label = "FPM transfer",
       p=3, c = cutoff, title = "Question 2")
