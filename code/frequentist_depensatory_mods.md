Simple depensatory models to all stocks
================

``` r
library(tidyverse)
```

    ## -- Attaching packages --------------------------------------- tidyverse 1.3.1 --

    ## v ggplot2 3.3.6     v purrr   0.3.4
    ## v tibble  3.1.7     v dplyr   1.0.7
    ## v tidyr   1.2.0     v stringr 1.4.0
    ## v readr   2.1.2     v forcats 0.5.1

    ## -- Conflicts ------------------------------------------ tidyverse_conflicts() --
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(ggh4x) #for some nicer facetting
library(here)
```

    ## here() starts at C:/Users/GLASERD/Desktop/projects/2022_depensation/salmon-depensation

#### Read in data

``` r
SR_data <- read.csv(here("data/salmon_productivity_compilation_jun2022.csv")) %>%
  select(-X, -stock.id) %>%
  mutate(species_stock = paste(species, stock)) #to help with looping later
```

## Fit Shepard models to all populaitons

Need to get rid of some populations that won’t fit.  
I think these aren’t fitting because of high recruitment variability,
which makes my simple init for the beta parm
(i.e. `max(sub_stock$spawners)/1.5`) break the model.

``` r
SR_data_filter <- filter(SR_data, !(species_stock %in% c("Chinook Andrew Creek", "Chinook Alsek/Klukshu",
                                                "Chinook Alsek-Klukshu",
                                                "Chinook East Fork South Fork Salmon River",
                                                "Chinook Entiat River", "Chinook Minam River",
                                                "Chinook Coweeman River"))) %>%
  mutate(spawners = ifelse(spawners==0|recruits==0, spawners+1, spawners),
         recruits = ifelse(spawners==0|recruits==0, recruits+1, recruits))
```

Then we can run the model

``` r
#build likliehood function
getNegLogLike <- function(parms){
  alpha <- parms[1]
  beta <- exp(parms[2])
  delta <- parms[3]
  sd <- exp(parms[4])
  
  predicted <- (alpha*sub_stock$spawners^delta)/(1+(beta*sub_stock$spawners^delta))
  NegLogLike <- -sum(log(dnorm(sub_stock$recruits, predicted, sd = sd)))
  return(NegLogLike)
}

#declare parms all pops will share
alpha <- 1
delta <- 1

shepard_results <- NULL #empty df to populate with results

for(i in unique(SR_data_filter$species_stock)){
  sub_stock <- filter(SR_data_filter, species_stock == i) %>%
      na.omit()
    
    #declare stock specific parms for beta and sd based on some heuristics 
    log_beta <- log(max(sub_stock$spawners)/1.5)
    log_sd <- log(sd(sub_stock$recruits/1.5)) 
    
    parms <- c(alpha, log_beta, delta, log_sd)
    
    #fit it
    fit1 <- optim(parms, getNegLogLike, method="Nelder-Mead", hessian=TRUE)
    fit2 <- optim(fit1$par, getNegLogLike, method="Nelder-Mead", hessian=TRUE)
    fit3 <- optim(fit2$par, getNegLogLike, method="Nelder-Mead", hessian=TRUE)
    

    shepard_results <- rbind(data.frame(i, unique(sub_stock$species), unique(sub_stock$stock),
                                        fit3$par[1], exp(fit3$par[2]), fit3$par[3], 
                                        exp(fit3$par[4])),
                             shepard_results)
}
```

And then clean up the table

``` r
colnames(shepard_results) <- c("species_stock", "species", "stock", "alpha", "beta", "delta", "sd")

shepard_results <- shepard_results %>% 
  mutate_at(4:7,funs(round(., 2)))
```

    ## Warning: `funs()` was deprecated in dplyr 0.8.0.
    ## Please use a list of either functions or lambdas: 
    ## 
    ##   # Simple named list: 
    ##   list(mean = mean, median = median)
    ## 
    ##   # Auto named with `tibble::lst()`: 
    ##   tibble::lst(mean, median)
    ## 
    ##   # Using lambdas
    ##   list(~ mean(., trim = .2), ~ median(., na.rm = TRUE))
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was generated.

Some of these `beta` parameters are tiny. They must have had a hard time
fitting. Let’s look at the data to try and get an idea of whats going
on. We’ll also include the pops that wouldn’t fit that we filtered out
earlier.

``` r
bad_stocks <- c(filter(shepard_results, beta<1) %>% pull(species_stock), 
                "Chinook Andrew Creek", "Chinook Alsek/Klukshu",
                                                "Chinook Alsek-Klukshu",
                                                "Chinook East Fork South Fork Salmon River",
                                                "Chinook Entiat River", "Chinook Minam River",
                                                "Chinook Coweeman River")


bad_stocks_data <- filter(SR_data, species_stock %in% bad_stocks)

max_facets <- 16

pages <- ceiling(length(bad_stocks)/max_facets)

for(i in 1:pages){
  sub_data <- filter(bad_stocks_data, species_stock %in% 
                               unique(bad_stocks_data$species_stock)[(((i-1)*max_facets)+1):(max_facets*i)])
  p <- ggplot(sub_data, aes(spawners, recruits)) +
    geom_point() +
        ggh4x::facet_wrap2(vars(species_stock), nrow=4, ncol=4, trim_blank=FALSE, 
                           scales="free")+
    labs(title = paste("Bad stocks", i, "of", pages))
  print(p)
}
```

![](frequentist_depensatory_mods_files/figure-gfm/check%20bad%20stocks-1.png)<!-- -->![](frequentist_depensatory_mods_files/figure-gfm/check%20bad%20stocks-2.png)<!-- -->![](frequentist_depensatory_mods_files/figure-gfm/check%20bad%20stocks-3.png)<!-- -->![](frequentist_depensatory_mods_files/figure-gfm/check%20bad%20stocks-4.png)<!-- -->

then let’s take a look at the delta estimates by species, for fits we
trust

``` r
good_fits <- filter(shepard_results, !(species_stock %in% bad_stocks))

ggplot(good_fits, aes(delta)) +
  geom_histogram() +
  facet_wrap(~species) +
  labs(title = "Histogram of delta parameter estimates form a depensatory beverton-holt")
```

    ## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.

![](frequentist_depensatory_mods_files/figure-gfm/compare%20deltas-1.png)<!-- -->
