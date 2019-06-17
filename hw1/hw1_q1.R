library(dplyr)

indicator <- function(x){
  I = ifelse(is.na(x), 1, 0) 
  return(I)
}

# read original data
data <- read.csv("wormwardata.csv")
# use only first visit in 1998 data
data1 <- filter(data, visit=="Visit 1, 1998")
data1 <- mutate(data1, is_na=indicator(totpar98))
# data with treated group in 1998
treated_data <- filter(data1, wgrp==1)

# create some vectors
schid <- data1$schid
unique_schid <- unique(schid)

# compute test statics with treated observations
is_na_treat <- treated_data$is_na
# obtain test statics
t_obs <- sum(is_na_treat)

# compute the distribution of test statics
# with Monte-Carlo approximation
n_sim = 1000
n <- length(data1$X)
is_na = data1$is_na
t_sim <- rep(NA, n_sim)
# iterate
for (b in 1:n_sim){
  treat_sch <- sample(unique_schid, size=25)
  treat_sim <- rep(NA, n)
  for (i in 1:n){
    treat_sim[i] <- ifelse(is.element(schid[i], treat_sch), 1, 0)
  }
  # compute test statics
  t_sim[b] <- sum(treat_sim * is_na)
}

# obtain p-value
pval <- mean(t_sim > t_obs)

# plot
hist(t_sim, xlab="Test Statics", main = paste("pval =", round(pval, 4)))
abline(v = t_obs, col='red', lwd = 1.5)
