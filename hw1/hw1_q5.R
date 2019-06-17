library(dplyr)

# indicator function
indicator <- function(x){
  I <- ifelse(x==1, 1, 0)
  return(I)
}

# read data
data <- read.csv("wormwardata.csv")
data1 <- filter(data, visit=="Visit 1, 1998")
data1 <- na.omit(data1)
data1 <- mutate(data1, treat=indicator(wgrp))

# treatment group data
data_treat <- filter(data1, treat == 1)
# control group data
data_contr <- filter(data1, treat == 0)

N <- length(data1$treat)
J <- length(unique(data1$schid))
J1 <- length(unique(data_treat$schid))
J2 <- length(unique(data_contr$schid))

# compute SATE
Y1b <- sum(data_treat$totpar98) * J / (N * J1)
Y0b <- sum(data_contr$totpar98) * J / (N * J0)
tau_hat <- Y1b - Y0b

# compute variance of treatment group
unique_treat_sch <- unique(data_treat$schid)
wY1 <- rep(NA, J1)
for (i in 1:J1){
  data_sch <- filter(data_treat, schid == unique_treat_sch[i])
  nj <- length(data_sch$totpar98)
  wY1[i] <- J * nj / N * mean(data_sch$totpar98)
}
VwY1 <- sum((wY1 - mean(wY1))^2) / (J1 - 1)

# compute variance of control group
unique_contr_sch <- unique(data_contr$schid)
wY0 <- rep(NA, J0)
for (i in 1:J0){
  data_sch <- filter(data_contr, schid == unique_contr_sch[i])
  nj <- length(data_sch$totpar98)
  wY0[i] <- J * nj / N * mean(data_sch$totpar98)
}
VwY0 <- sum((wY0 - mean(wY0))^2) / (J0 - 1)

# variance of SATE
Vtau <- VwY1 / J1 + VwY0 / J0
# standard deviation of SATE
sdtau <- sqrt(Vtau)

# p-value
p1 <- pnorm(tau_hat, mean = 0, sd = sdtau)
pval <- 2 * min(p1, 1 - p1)

# asymptotic 95% confidence interval
alpha = 0.05
ub <- tau_hat - qnorm(alpha/2) * sdtau
lb <- tau_hat + qnorm(alpha/2) * sdtau
