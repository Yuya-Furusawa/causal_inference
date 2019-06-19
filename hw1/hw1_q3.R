library(dplyr)

# indicator function
indicator <- function(x){
  I <- ifelse(x==1, 1, 0)
  return(I)
}

# fisher's confidence interaval : difference-in-means
fisher_CI_dim <- function(Y, treat, schid, n_sim, taus){
  pval <- rep(NA, length(taus))
  unique_schid <- unique(schid)
  for (i in 1:length(taus)){
    
    # create potential outcome
    Y01 <- cbind(Y,Y)
    Y01[treat==1,1] <- Y[treat==1] - taus[i]
    Y01[treat==0,2] <- Y[treat==0] + taus[i]
    
    # conduct a test
    t_obs <- abs(mean(Y01[treat==1,2]) - mean(Y01[treat==0,1]))
    t_sim <- rep(NA, n_sim)
    for (b in 1:n_sim){
      treat_sch <- sample(unique_schid, size=25)
      treat_sim <- rep(NA, length(Y))
      treat_sim <- if_else(is.element(schid, treat_sch), 1, 0)
      t_sim[b] <- abs(mean(Y01[treat_sim==1,2]) - mean(Y01[treat_sim==0,1]))
    }
    
    # obtain p-value
    pval[i] <- mean(t_sim > t_obs)
  }
  return(pval)
}

# fisher's confidence interaval : rank-sum
fisher_CI_rank <- function(Y, treat, schid, n_sim, taus){
  pval <- rep(NA, length(taus))
  unique_schid <- unique(schid)
  for (i in 1:length(taus)){
    
    # create potential outcome
    Y01 <- cbind(Y,Y)
    Y01[treat==1,1] <- Y[treat==1] - taus[i]
    Y01[treat==0,2] <- Y[treat==0] + taus[i]
    
    # observed statistics
    y_rank <- rank(Y, ties.method = "random") - (length(Y)+1)/2
    t_obs <- sum(treat * y_rank)
    
    # conduct a test
    t_sim <- rep(NA, n_sim)
    for (b in 1:n_sim){
      treat_sch <- sample(unique_schid, size=25)
      treat_sim <- rep(NA, length(Y))
      treat_sim <- if_else(is.element(schid, treat_sch), 1, 0)
      Y_sim <- rep(NA, length(Y))
      Y_sim <- if_else(treat_sim==1, Y01[,2], Y01[,1])
      y_rank_sim <- rank(Y_sim, ties.method = "random") - (length(Y)+1)/2
      t_sim[b] <- sum(treat_sim * y_rank_sim)
    }
    
    # obtain p-value
    pval[i] <- mean(t_sim > t_obs)
  }
  return(pval)
}

# read data
data <- read.csv("wormwardata.csv")
data1 <- filter(data, visit=="Visit 1, 1998")
data1 <- na.omit(data1)
data1 <- mutate(data1, treat=indicator(wgrp))

set.seed(1234)
taus <- seq(0, 0.25, by=0.005)
fr_ci_dim <- fisher_CI_dim(data1$totpar98, data1$treat, data1$schid, 500, taus)
fr_ci_rank <- fisher_CI_rank(data1$totpar98, data1$treat, data1$schid, 500, taus)

alpha <- 0.05
accept_dim <- fr_ci_dim >= alpha/2 & fr_ci_dim <= (1 - alpha/2)
CI_dim <- taus[accept_dim]
accept_rank <- fr_ci_rank >= alpha/2 & fr_ci_rank <= (1 - alpha/2)
CI_rank <- taus[accept_rank]

# plot
par(mfrow = c(2,1))
plot(taus, fr_ci_dim, xlab="taus", ylab="pval(tau)", pch=21, bg=ifelse(accept_dim==1, 'black', 'gray'))
abline(h = 1 - alpha/2, col = 'red')
abline(h = alpha/2, col = 'red')

plot(taus, fr_ci_rank, xlab="taus", ylab="pval(tau)", pch=21, bg=ifelse(accept_rank==1, 'black', 'gray'))
abline(h = 1 - alpha/2, col = 'red')
abline(h = alpha/2, col = 'red')
