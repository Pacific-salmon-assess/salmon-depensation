#simple shepard funs fit to all of the pops in Dan's data

library(tidyverse)

SR_data <- read.csv("data/salmon_productivity_compilation_jun2022.csv") %>%
  select(-X, -stock.id)

##Old way of coding a nonlinear from cahill days##
SR_test <- filter(SR_data, stock == "Atnarko") #subset 1 stock and try to solve

#build NLL function - Same one as from Myers 1995
#start parms
alpha <- 1
log_beta <- log(max(SR_test$spawners)/1.5)
delta <- 1
log_sd <- log(sd(SR_test$spawners/1.5)) 

parms <- c(alpha, log_beta, delta, log_sd)

getNegLogLike <- function(parms){
  alpha <- parms[1]
  beta <- exp(parms[2])
  delta <- parms[3]
  sd <- exp(parms[4])
  
  predicted <- (alpha*SR_test$spawners^delta)/(1+(beta*SR_test$spawners^delta))
  NegLogLike <- -sum(log(dnorm(SR_test$recruits, predicted, sd = sd)))
  return(NegLogLike)
}

#fit it several times. Not really sure this is "good" because it might just hone in
  #on a local minima more! lol 
fit1 <- optim(parms, getNegLogLike, method="Nelder-Mead", hessian=TRUE)
fit2 <- optim(fit1$par, getNegLogLike, method="Nelder-Mead", hessian=TRUE)
fit3 <- optim(fit2$par, getNegLogLike, method="Nelder-Mead", hessian=TRUE)

fit3$par

#ok proof of concept is done. time to loop it for everyon

if(FALSE){
##try another way with nls() wow wtf I cant figure this out 
#use linearized ricker to get a reasonable starting alpha - is this apples to apples though?
summary(lm(log(recruits/spawners) ~ spawners, data = SR_test))

fit <- nls(recruits ~ (alpha*spawners*exp(delta))/(1+(beta*spawners*exp(delta))),
           data = SR_test,
           start = list(alpha = 0.4, beta = max(SR_test$spawners)/1.5, delta = 1), 
           algorithm = "plinear")
}