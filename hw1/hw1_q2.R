library(dplyr)

# indicator function
indicator <- function(x){
  I <- ifelse(x==1, 1, 0)
  return(I)
}

# fisher's randomize test : difference-in-means
fisher_randomize_dim <- function(Y, treat, schid, n_sim){

  # observed statistics
  t_obs <- abs(mean(Y[treat==1]) - mean(Y[treat==0]))
  
  # Monte-Carlo approximation
  t_sim <- rep(NA, n_sim)
  for (b in 1:n_sim){
    unique_schid <- unique(schid)
    treat_sch <- sample(unique_schid, size=25)
    treat_sim <- rep(NA, length(Y))
    for (i in 1:length(Y)){
      treat_sim[i] <- ifelse(is.element(schid[i], treat_sch), 1, 0)
    }
    t_sim[b] <- abs(mean(Y[treat_sim==1]) - mean(Y[treat_sim==0]))
  }
  
  # obtain p-value
  pval <- mean(t_sim > t_obs)
  
  return(list(pval=pval, t_obs=t_obs, t_sim=t_sim))
}

# fisher's randomization test : rank-sum
fisher_randomize_rank <- function(Y, treat, schid, n_sim) {
  
  # observed statistics
  y_rank <- rank(Y, ties.method = "random") - (length(Y)+1)/2
  t_obs <- sum(treat * y_rank)
  
  # Monte-Carlo approximation
  t_sim <- rep(NA, n_sim)
  for (b in 1:n_sim){
    unique_schid <- unique(schid)
    treat_sch <- sample(unique_schid, size=25)
    treat_sim <- rep(NA, length(Y))
    for (i in 1:length(Y)){
      treat_sim[i] <- ifelse(is.element(schid[i], treat_sch), 1, 0)
    }
    t_sim[b] <- sum(treat_sim * y_rank)
  }
  
  # obtain p-value
  pval <- mean(t_sim > t_obs)
  
  return(list(pval=pval, t_obs=t_obs, t_sim=t_sim))
}

# read data
data <- read.csv("wormwardata.csv")
data1 <- filter(data, visit=="Visit 1, 1998")
data1 <- na.omit(data1)
data1 <- mutate(data1, treat=indicator(wgrp))

set.seed(1234)
fr1 <- fisher_randomize_dim(data1$totpar98, data1$treat, data1$schid, 1000)
fr2 <- fisher_randomize_rank(data1$totpar98, data1$treat, data1$schid, 1000)

# plot
par(mfrow = c(1,2))
hist(fr1$t_sim, xlab="difference-in-means", main=paste("pval =", round(fr1$pval, 4)))
abline(v = fr1$t_obs, col = 'red', lwd = 1.5)
hist(fr2$t_sim, xlab="normalized rank", main=paste("pval =", round(fr2$pval, 4)))
abline(v = fr2$t_obs, col = 'red', lwd = 1.5)
