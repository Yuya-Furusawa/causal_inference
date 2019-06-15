fisher_CI <- function(Y, treat, n_sim, taus){
  pval <- rep(NA, length(taus))
  for (i in 1:length(taus)){
    
    # create potential outcome
    Y01 <- cbind(Y,Y)
    Y01[treat==1,1] <- Y[treat==1] - taus[i]
    Y01[treat==0,2] <- Y[treat==0] + taus[i]
    
    # conduct a test
    t_obs <- abs(mean(Y01[treat==1,2]) - mean(Y01[treat==0,1]))
    t_sim <- rep(NA, n_sim)
    for (b in 1:n_sim){
      treat_sim <- sample(treat, size = length(treat))
      t_sim[b] <- abs(mean(Y01[treat_sim==1,2]) - mean(Y01[treat_sim==0,1]))
    }
    
    # obtain p-value
    pval[i] <- mean(t_sim > t_obs)
  }
  return(pval)
}