library(dplyr)

data <- readRDS("enos2016.rds")
data <- data %>%
  mutate(treated=ifelse(demo.distance <= 500, 1, 0))
n <- length(data$id)

# Question1-c
data_q1_treated <- data %>%
  filter(treated == 1)
data_q1_control <- data %>%
  filter(treated == 0)
treated = c(mean(data_q1_treated$vote1996),
            mean(data_q1_treated$vote2000),
            mean(data_q1_treated$vote2004))
control = c(mean(data_q1_control$vote1996),
            mean(data_q1_control$vote2000),
            mean(data_q1_control$vote2004))
plot(c(1996, 2000, 2004), treated, ylim = c(0.4,0.8),
     type='l', xlab="year", ylab="mean prob. turnout", col="red", main="prallel trend")
par(new=TRUE)
plot(c(1996, 2000, 2004), control, ylim = c(0.4,0.8),
     type='l', xlab="", ylab="", col="black")
legend("bottomright", legend = c("treated", "control"), col = c("red", "black"), lwd=c(1,1))

# Question1-d
delta_1 <- data_q1_treated$vote2000 - data_q1_treated$vote1996
delta_0 <- data_q1_control$vote2000 - data_q1_control$vote1996
t.test(delta_1, delta_0, alternative="two.sided", var.equal=T)

# Question1-e
time <- c(rep(0,n), rep(1,n))
z <- c(data$treated, data$treated)
y <- c(data$vote2000, data$vote2004)
model.did = lm(y ~ time + z)
summary(model.did)


# Question2-b
##(1)
data_q2 <- data %>%
  mutate(white=ifelse(whitename >= 0.975, 1, 0)) %>%
  mutate(black=ifelse(blackname >= 0.975, 1, 0))
data_q2_white <- data_q2 %>% filter(white==1)
data_q2_black <- data_q2 %>% filter(black==1)
t_white <- data_q2_white %>% filter(treated==1)
c_white <- data_q2_white %>% filter(treated==0)
t_black <- data_q2_black %>% filter(treated==1)
c_black <- data_q2_black %>% filter(treated==0)
treated_white = c(mean(t_white$vote1996),
                  mean(t_white$vote2000),
                  mean(t_white$vote2004))
control_white = c(mean(c_white$vote1996),
                  mean(c_white$vote2000),
                  mean(c_white$vote2004))
treated_black = c(mean(t_black$vote1996),
                  mean(t_black$vote2000),
                  mean(t_black$vote2004))
control_black = c(mean(c_black$vote1996),
                  mean(c_black$vote2000),
                  mean(c_black$vote2004))
plot(c(1996, 2000, 2004), treated_white, ylim = c(0.3,0.9),
     type='l', xlab="year", ylab="mean prob. turnout", col="red", main="parallel trend : white")
par(new=TRUE)
plot(c(1996, 2000, 2004), control_white, ylim = c(0.3,0.9),
     type='l', xlab="", ylab="", col="black")
legend("bottomright", legend = c("treated", "control"), col = c("red", "black"), lwd=c(1,1))
par(new=False)
plot(c(1996, 2000, 2004), treated_black, ylim = c(0.3,0.9),
     type='l', xlab="year", ylab="mean prob. turnout", col="red", main="parallel trend : black")
par(new=TRUE)
plot(c(1996, 2000, 2004), control_black, ylim = c(0.3,0.9), type='l', xlab="", ylab="", col="black")
legend("bottomright", legend = c("treated", "control"), col = c("red", "black"), lwd=c(1,1))

##(2)
n_white <- length(data_q2_white$id)
time_white <- c(rep(0,n_white), rep(1,n_white))
z_white <- c(data_q2_white$treated, data_q2_white$treated)
y_white <- c(data_q2_white$vote2000, data_q2_white$vote2004)
model.did_white = lm(y_white ~ time_white + z_white)
summary(model.did_white)

n_black <- length(data_q2_black$id)
time_black <- c(rep(0,n_black), rep(1,n_black))
z_black <- c(data_q2_black$treated, data_q2_black$treated)
y_black <- c(data_q2_black$vote2000, data_q2_black$vote2004)
model.did_black = lm(y_black ~ time_black + z_black)
summary(model.did_black)

# Question2-(e)
data_q2 <- data_q2 %>%
  mutate(log.medianincome=log(medianincome)) %>%
  mutate(log.prior.ave.value=log(prior.avg.value)) %>%
  mutate(gender.dummy=ifelse(gender=="M", 1, 0))

pi <- sum(data_q2$treated) / n

logi <- glm(treated ~ pid + gender.dummy + age + log.medianincome + log.prior.ave.value + pctblack, data=data_q2, family=binomial)
ps <- logi$fitted.values

total = 0
for (i in 1:n){
  total = total + (data_q2$vote2004[i] - data_q2$vote2000[i]) / pi * (data_q2$treated[i] - ps[i]) / (1 - ps[i])
}
tau_ht = total/n

total = 0
for (i in 1:n){
  total = total + ((data_q2$vote2004[i] - data_q2$vote2000[i]) / pi * (data_q2$treated[i] - ps[i]) / (1 - ps[i]) - tau_ht)^2
}
tau_ht_sd = sqrt(total/n)

q2e <- matrix(0, 1, 2)
colnames(q2e) <- c("estimated tau", "standard error") 
q2e[1,] <- c(tau_ht, tau_ht_sd)
q2e


# Question3-(c)
n1 <- sum(data_q2$treated)
n0 <- n - n1
tau1 = sum(data_q2$treated * (data_q2$vote2004 - data_q2$vote2000)) / n1 - sum((rep(1,n) -data_q2$treated) * (data_q2$vote2004 - data_q2$vote2000)) / n0
tau2 = sum(data_q2$treated * (data_q2$vote2004 - data_q2$vote1996)) / n1 - sum((rep(1,n) -data_q2$treated) * (data_q2$vote2004 - data_q2$vote1996)) / n0
beta_hat = (tau1 + tau2) / 2

