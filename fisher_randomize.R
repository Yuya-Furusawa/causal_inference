fisher_randomize_dim <- function(Y, treat, n_sim) {
  
  # observed statistics
  t_obs <- abs(mean(Y[treat==1]) - mean(Y[treat==0]))
  
  # Monte-Carlo approximation
  t_sim <- rep(NA, n_sim)
  for (b in 1:n_sim){
    treat_sim <- sample(treat, size=length(treat))
    t_sim[b] <- abs(mean(Y[treat_sim==1]) - mean(Y[treat_sim==0]))
  }
  
  # obtain p-value
  pval <- mean(t_sim > t_obs)
  
  return(list(pval=pval, t_obs=t_obs, t_sim=t_sim))
}

fisher_randomize_rank <- function(Y, treat, n_sim) {
  
  # observed statistics
  y_rank <- rank(Y, ties.method = "random") - (length(Y)+1)/2
  t_obs <- sum(treat * y_rank)
  
  # Monte-Carlo approximation
  t_sim <- rep(NA, n_sim)
  for (b in 1:n_sim){
    treat_sim <- sample(treat, size=length(treat))
    t_sim[b] <- sum(treat_sim * y_rank)
  }
  
  # obtain p-value
  pval <- mean(t_sim > t_obs)
  
  return(list(pval=pval, t_obs=t_obs, t_sim=t_sim))
}