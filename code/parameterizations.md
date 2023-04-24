Various spawner recruit relationships for depensatory dynamics
================

# Background

Here we explore various parameterizations of stock recruit relationships
(SRRs) that show evidence of depensation. There are many
parameterizaitons for both Ricker and Beverton-Holt type spawner recruit
relationships (SRRs) but we will focus on some basics. Most of these
derivations arise from a need for how to treat the depensatory
parameter, with some having a goal of making the depensatory parameter
biologically relevant.

## Ricker based SRRs

There are 3 types we will explore here.  
- original [Ricker](https://cdnsciencepub.com/doi/10.1139/f54-039) from
1954  
- depensatory Ricker ([Quinn & Deriso,
1999](https://books.google.ca/books/about/Quantitative_Fish_Dynamics.html?id=5FVBj8jnh6sC&redir_esc=y),
p. 99)  
- Saila-Lorda (as described in [Perälä & Kuparinen,
2017](http://dx.doi.org/10.1098/rspb.2017.1284))

To understand these models we’ll describe their structure, simulate data
based on the Saila-Lorda then fit each SRR to the simulated data.

#### Original Ricker

$$ R = \alpha Se^{-\beta S} $$  
where $\alpha$ is the slope near the origin and $\beta$ is the density
dependent parameter. There is no term for density dependence in this
equation, but we will use it to compare. The Ricker SRR is also popular
because it can be made linear:  
$$ log(R/S) = \alpha - \beta S $$

#### Quinn & Deriso’s 4 parameter depensatory Ricker SRR

$$ R = \alpha Se^{\beta S - \delta S^{\gamma }} $$ The authors say that
this model is depensatory when
$S^{\gamma -1} < \beta / (\gamma \delta)$, which isn’t very intuitive,
which brings us to the next model.

#### Salia-Lorda SRR

$$ R_{t} = k(\frac{S_{t}}{S_{k}})^{c}e^{c (1-S_{t}/S_{k})} $$ where
$R_{t}$ and $S_{t}$ are vectors of recruits and spawners at time *t*,
*k* is the maximum number of recruits produced when $S_{t} = S_{k}$ and
$S_{k}$ is is the spawning stock at which the number of recruits is *k*.
*c* is the depensatory parameter; *c* = 1 produces the classic Ricker,
*c* \> 1 shows depensatory dynamics, and *c* \< 1 shows
hypercompensaiton.

## Simulation

We’ll simulate some data from a Salia-Lorda to play with. FEEL FREE to
toggle the parameters (i.e. S, Sk, k, c, sigma), they will cascade
through the robust code

``` r
S <- runif(100, 10, 1000) #vector of random spawners
Sk <- 600 
k <- 200
c <- 3
sigma <- rnorm(length(S), 0, 50) #random white noise error
```

![](parameterizations_files/figure-gfm/write%20fun%20plot%20sim-1.png)<!-- -->

## Fit simulated data

Can fit a Ricker here the “old fashioned” way by making it linear.

    ##                 Estimate   Std. Error   t value     Pr(>|t|)
    ## (Intercept) -0.492495722 0.1045571842 -4.710300 9.058274e-06
    ## spawners    -0.001349195 0.0001626803 -8.293539 1.082784e-12

now we can use optim to fit a non linear Quinn & Deriso model  
<table>
<tbody>
<tr>
<td style="text-align:left;">
alpha
</td>
<td style="text-align:left;">
0.6111
</td>
</tr>
<tr>
<td style="text-align:left;">
beta
</td>
<td style="text-align:left;">
-0.0013
</td>
</tr>
<tr>
<td style="text-align:left;">
delta
</td>
<td style="text-align:left;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
gamma
</td>
<td style="text-align:left;">
1
</td>
</tr>
<tr>
<td style="text-align:left;">
sd
</td>
<td style="text-align:left;">
2.3026
</td>
</tr>
</tbody>
</table>
Now we’ll refit the Salia-Lorda and plot the results
<table>
<tbody>
<tr>
<td style="text-align:left;">
Sk
</td>
<td style="text-align:left;">
617.97
</td>
</tr>
<tr>
<td style="text-align:left;">
k
</td>
<td style="text-align:left;">
190.13
</td>
</tr>
<tr>
<td style="text-align:left;">
c
</td>
<td style="text-align:left;">
2.04
</td>
</tr>
<tr>
<td style="text-align:left;">
sd
</td>
<td style="text-align:left;">
3.93
</td>
</tr>
</tbody>
</table>

Now we’ll generate model predictions from the parameter fits

Then plot it
![](parameterizations_files/figure-gfm/plot%20fits-1.png)<!-- -->

We could add variations of the Beverton-Holt, which are generally easier
to fit. These could include:  
- classic Beverton-Holt without depensation  
- Myers’ Beverton-Holt from the ’95
[paper](https://www.science.org/doi/10.1126/science.269.5227.1106) in
*Science*  
- The Beverton-Holt fit in the Perala et al
[paper](https://royalsocietypublishing.org/doi/10.1098/rsbl.2021.0439)  
- [Lierman &
Hilborn’s](https://cdnsciencepub.com/doi/abs/10.1139/f97-105) stepwise
way of fitting a depensatory Beverton-Holt with meaningful parameters
