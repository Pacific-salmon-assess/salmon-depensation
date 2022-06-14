Various stock recruit relationships for depensatory dynamics
================

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ purrr   0.3.4
    ## ✔ tibble  3.1.7     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

## Background

Here we explore various paramaterizations of stock recruit relationships
(SRRs) that show evidence of depensation. There are many
parameterizaitons for both Ricker and Beverton-Holt type SRRs but we
will focus on some basics. Most of these derivations arise from a need
for how to treat the depensatory parameter, with some having a goal of
making the depensatory parameter “biologically relevant”.

### Ricker based SRRs

There are 3 types we will explore here. - original
[Ricker](https://cdnsciencepub.com/doi/10.1139/f54-039) from 1954  
- depensatory Ricker ([Quinn & Deriso,
1999](https://books.google.ca/books/about/Quantitative_Fish_Dynamics.html?id=5FVBj8jnh6sC&redir_esc=y),
p. 99) - Saila-Lorda (as described in [Perala & Kuparinen,
2017](http://dx.doi.org/10.1098/rspb.2017.1284))

To understand these models we’ll describe their structure, simulate data
based on the Saila-Lorda then fit each SRR to the simulated data.

#### Original Ricker

![ R = \\alpha Se^{-\\beta S} ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20R%20%3D%20%5Calpha%20Se%5E%7B-%5Cbeta%20S%7D%20 " R = \alpha Se^{-\beta S} ")

  
where
![\\alpha](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Calpha "\alpha")
is the slope near the origin and
![\\beta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta "\beta")
is the density dependent parameter. There is no term for density
dependence in this equation, but we will use it to compare. The Ricker
SRR is also popular because it can be transformed to be linear:  

![ log(R/S) = \\alpha - \\beta S ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20log%28R%2FS%29%20%3D%20%5Calpha%20-%20%5Cbeta%20S%20 " log(R/S) = \alpha - \beta S ")

#### Quinn & Deriso’s 4 parameter depensatory Ricker SRR

![ R = \\alpha Se^{\\beta S - \\delta S^{\\gamma }} ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20R%20%3D%20%5Calpha%20Se%5E%7B%5Cbeta%20S%20-%20%5Cdelta%20S%5E%7B%5Cgamma%20%7D%7D%20 " R = \alpha Se^{\beta S - \delta S^{\gamma }} ")

The authors say that this model is depensatory when
![S^{\\gamma -1} \< \\beta / (\\gamma \\delta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S%5E%7B%5Cgamma%20-1%7D%20%3C%20%5Cbeta%20%2F%20%28%5Cgamma%20%5Cdelta%29 "S^{\gamma -1} < \beta / (\gamma \delta)"),
which isn’t very intuitive, which brings us to the next model.

#### Salia-Lorda SRR

![ R\_{t} = k(\\frac{S\_{t}}{S\_{k}})^{C}e^{C (1-S\_{t}/S\_{k})} ](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%20R_%7Bt%7D%20%3D%20k%28%5Cfrac%7BS_%7Bt%7D%7D%7BS_%7Bk%7D%7D%29%5E%7BC%7De%5E%7BC%20%281-S_%7Bt%7D%2FS_%7Bk%7D%29%7D%20 " R_{t} = k(\frac{S_{t}}{S_{k}})^{C}e^{C (1-S_{t}/S_{k})} ")

where
![R\_{t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_%7Bt%7D "R_{t}")
and
![S\_{t}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_%7Bt%7D "S_{t}")
are vectors of recruits and spawners at time *t*, *k* is the maximum
number of recruits produced when
![S\_{t} = S\_{k}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_%7Bt%7D%20%3D%20S_%7Bk%7D "S_{t} = S_{k}")
and
![S\_{k}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;S_%7Bk%7D "S_{k}")
is is the spawning stock at which the number of recruits is *k*. *c* is
the depensatory parameter; *c* = 1 produces the classic Ricker, *c* \> 1
shows depensatory dynamics, and *c* \< 1 shows hypercompensaiton.

## Simulation

We’ll simulate some data to play with. FEEL FREE to toggle the
parameters (i.e. S, Sk, k, c, sigma), they will cascade through the
robust code

``` r
sal_lor <- function(S, Sk, k, c, sigma){
  exp_rec <- ((k*(S/Sk)^c)*exp(c*(1-S/Sk))+sigma)
  return(exp_rec)
}
  
S <- runif(100, 0, 1000)
Sk <- 600
k <- 200
c <- 3
sigma <- rnorm(length(S), 0, 10)

exp_rec <- sal_lor(S, Sk, k, c, sigma)

plot(S, exp_rec)
```

![](parameterizations_files/figure-gfm/write%20fun-1.png)<!-- -->

``` r
sim_stock <- data.frame(stock = S, recruits = exp_rec)
```

## Fit simulated data

Can fit a ricker here the “old fashinoed way by making it linear.

``` r
rick_fit <- lm(log(recruits/stock)~stock, data = sim_stock)
```

    ## Warning in log(recruits/stock): NaNs produced

``` r
summary(rick_fit)$coefficients
```

    ##                 Estimate   Std. Error   t value     Pr(>|t|)
    ## (Intercept) -0.929777091 0.1149746440 -8.086801 1.975182e-12
    ## stock       -0.000736401 0.0001969179 -3.739636 3.150275e-04

now we can use optim to fit a non linear Quinn & Deriso model

``` r
QDLogLike <- function(parms){
  alpha <- parms[1]
  beta <- parms[2]
  delta <- parms[3]
  gamma <- parms[4]
  sd <- exp(parms[5])
  
  predicted <- alpha*sim_stock$stock*exp((beta*sim_stock$stock)-(delta*sim_stock$stock^gamma))
  
  NegLogLike <- -sum(log(dnorm(sim_stock$recruits, predicted, sd= sd)))
  return(NegLogLike)
}

#declare start parms - guess based on previous fit 
alpha <- exp(summary(rick_fit)$coefficients[1])
beta <- -summary(rick_fit)$coefficients[2]
delta <- 1
gamma <- 1
sd <- log(10)

parms <- list(alpha, beta, delta, gamma, sd)

QD_fit <- optim(parms, QDLogLike, method="Nelder-Mead", hessian=TRUE)

QD_fit$par
```

    ## [1]  0.8406881 -0.3514120 -0.7352116  0.2919089  4.9556950

Now we’ll refit the Salia-Lorda

``` r
SLLogLike <- function(parms){
  Sk <- parms[1]
  k <- parms[2]
  c <- parms[3]
  sd <- exp(parms[4])
  
  predicted <- ((k*(sim_stock$stock/Sk)^c)*exp(c*(1-sim_stock$stock/Sk)))  
  NegLogLike <- -sum(log(dnorm(sim_stock$recruits, predicted, sd= sd)))
  return(NegLogLike)
}
#Sk <- 1000 #tester - this is finnecky! and changes results - must be hitting local minima
parms <- list(Sk, k, c, sd) #can recycle from above

SL_fit <- optim(parms, SLLogLike, method="Nelder-Mead", hessian=TRUE)
cbind(c("SK", "k", "c", "sd"), SL_fit$par)
```

    ##      [,1] [,2]              
    ## [1,] "SK" "596.168642310362"
    ## [2,] "k"  "198.965534625435"
    ## [3,] "c"  "2.82656189534481"
    ## [4,] "sd" "2.24743231071727"

Now we’ll gnerate model predictions from the parameter fits

``` r
r_alpha <- exp(summary(rick_fit)$coefficients[1])
r_beta <- summary(rick_fit)$coefficients[2]

rick_preds <- r_alpha*sim_stock$stock*exp(r_beta*sim_stock$stock)
rick_preds <- data.frame(fit = rick_preds, S = S) %>%
  arrange(S)

QD_alpha <- QD_fit$par[1]
QD_beta <- QD_fit$par[2]
QD_delta <- QD_fit$par[3]
QD_gamma <- QD_fit$par[4]

QD_preds <- QD_alpha*sim_stock$stock*exp((QD_beta*sim_stock$stock)-(QD_delta*sim_stock$stock^QD_gamma))
QD_preds <- data.frame(fit = QD_preds, S = S) %>%
  arrange(S) 
SL_Sk <- SL_fit$par[1]
SL_k <- SL_fit$par[2]
SL_c <- SL_fit$par[3]

SL_preds <- ((SL_k*(sim_stock$stock/Sk)^c)*exp(c*(1-sim_stock$stock/Sk)))  

SL_preds <- data.frame(fit = SL_preds, S = S) %>%
  arrange(S)
```

Then plot it

``` r
plot(S, exp_rec)
lines(fit~S, data = rick_preds, col = "blue")
lines(fit~S, data = SL_preds, col = "green")
lines(fit~S, data = QD_preds, col = "red") #Something's messed up with the QD model!
```

![](parameterizations_files/figure-gfm/plot%20fits-1.png)<!-- -->
