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

# RD plot : narrow
rdplot(y=stan_data$narrow, x=stan_data$pop,
       x.label = "population", y.label = "irregurarity",
       title="narrow", p=3)

# RD plot : broad
rdplot(y=stan_data$broad, x=stan_data$pop,
       x.label = "population", y.label = "irregurarity",
       title="broad", p=3)

# Compute Density
n <- length(stan_data$pop)
b <- 2 * sd(stan_data$pop) / sqrt(n)
Xdisc <- floor(stan_data$pop/b) * b + b/2
m <- (max(Xdisc) - min(Xdisc)) / min(diff(sort(unique(Xdisc)))) + 1
Xj <- min(Xdisc) + (0:m * min(diff(sort(unique(Xdisc)))))
Yj <- table(factor(round(Xdisc,3), levels=round(Xj,3))) / (n*b)
rdplot(y=Yj, x=Xj, x.label = "standardized variable", y.label = "Normalized count")

# RD Regression
fit <- rdrobust(y=stan_data$broad, x=stan_data$pop, p=1, kernel="tri")
res_tab <- cbind(fit$coef, fit$se)
rownames(res_tab) <- rownames(fit$coef)
colnames(res_tab) <- c("estimate", "se")
print(res_tab)

fit <- rdrobust(y=stan_data$narrow, x=stan_data$pop, p=1, kernel="tri")
res_tab <- cbind(fit$coef, fit$se)
rownames(res_tab) <- rownames(fit$coef)
colnames(res_tab) <- c("estimate", "se")
print(res_tab)
