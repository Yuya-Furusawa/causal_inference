library(dplyr)

# read data
data <- read.csv("corruption.csv")

# create standardized data
cutoffs <- c(10188, 13584, 16980, 23772, 30564, 37356, 44148)
cdata <- vector("list", length(cutoffs))
for (i in 1:length(cutoffs)){
  cdata[[i]] <- filter(data, pop >= cutoffs[i] - 1000, pop <= cutoffs[i] + 1000)
  cdata[[i]]$pop <- cdata[[i]]$pop - cutoffs[i]
}
stan_data <- rbind(cdata[[1]], cdata[[2]])
for (i in 3:length(cutoffs)){
  stan_data <- rbind(stan_data, cdata[[i]])
}

fit_frd <- rdrobust(y=stan_data$broad, x=stan_data$pop, fuzzy=stan_data$fpm, p=1, kernel="tri")
summary(fit_frd)

fit_frd <- rdrobust(y=stan_data$narrow, x=stan_data$pop, fuzzy=stan_data$fpm, p=1, kernel="tri")
summary(fit_frd)

# 2SLS
indicator <- function(x){
  I <- ifelse(x >= 0, 1, 0)
  return(I)
}

stan_data <- mutate(stan_data, above=indicator(pop))

result <- lm(stan_data$fpm ~ stan_data$pop * stan_data$above)
D_hat <- stan_data$fpm - result1$residuals  # fitted value of D

result_broad <- lm(stan_data$broad ~ D_hat * stan_data$pop)
result_narrow <- lm(stan_data$narrow ~ D_hat * stan_data$pop)