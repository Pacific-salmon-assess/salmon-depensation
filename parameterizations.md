Various stock recruit relationships for depensatory dynamics
================

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

We’ll simulate some data to play with.

``` r
sal_lor <- function(S, Sk, k, c, sigma){
  exp_rec = ((k*(S/Sk)^c)*exp(c*(1-S/Sk))+sigma)
  return(exp_rec)
}
  
S <- runif(100, 0, 1000)
Sk <- 500
k <- 100
c <- 3
sigma <- rnorm(length(S), 0, 10)

exp_rec <- sal_lor(S, Sk, k, c, sigma)

plot(S, exp_rec)
```

![](parameterizations_files/figure-gfm/write%20fun-1.png)<!-- -->
