#simple shepard funs fit to all of the pops in Dan's data

library(tidyverse)

SR_data <- read.csv("data/salmon_productivity_compilation_jun2022.csv") %>%
  select(-X, -stock.id)

#build NLL function - Same one as from Myers 1995
#start parms
alpha <- 1
beta <- log(30000)
delta <- 1
#sd <- log(20)

parms <- c(alpha, beta, delta)

getNegLogLike <- function(parms){
  alpha <- parms[1]
  beta <- exp(parms[2])
  delta <- parms[3]
  #sd <- exp(parms[4])
  
  predicted <- (alpha*SR_test$spawners^delta)/(1+(beta*SR_test$spawners^delta))
  NegLogLike <- -sum(log(dnorm(observed, predicted)))
  return(NegLogLike)
}

SR_test <- filter(SR_data, stock == "Atnarko") #subset 1 stock and try to solve

predicted <- (alpha*SR_test$spawners^delta)/(1+(SR_test$spawners^delta/exp(beta)))
observed  <- SR_test$recruits

optim(parms, getNegLogLike, method="Nelder-Mead", hessian=TRUE)

#i dont think the code from this old exercise is right for this application. I should do 
#this in nlme or something. 