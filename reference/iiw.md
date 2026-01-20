# Given a proportional hazards model for visit intensities, compute inverse-intensity weights.

For a longitudinal dataset subject to irregular observation, use a Cox
proportional hazards model for visit intensities to compute inverse
intensity weights

## Usage

``` r
iiw(phfit, data, id, time, first)
```

## Arguments

- phfit:

  coxph object for the visit process

- data:

  The dataset featuring longitudinal data subject to irregular
  observation for which inverse-intensity weights are desired

- id:

  character string indicating which column of the data identifies
  subjects

- time:

  character string indicating which column of the data contains the time
  at which the visit occurred

- first:

  logical variable. If TRUE, the first observation for each individual
  is assigned an intensity of 1. This is appropriate if the first visit
  is a baseline visit at which recruitment to the study occurred; in
  this case the baseline visit is observed with probability 1.

## Value

A vector of inverse-intensity weights for each row of the dataset. The
first observation for each subject is assumed to have an intensity of 1.

## See also

Other iiw:
[`iiw.weights()`](https://epullenayegum.github.io/IrregLong/reference/iiw.weights.md),
[`iiwgee()`](https://epullenayegum.github.io/IrregLong/reference/iiwgee.md)

## Examples

``` r
library(nlme)
library(survival)
library(geepack)
library(data.table)
data(Phenobarb)
Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
data <- Phenobarb
data <- data[data$event==1,]
data$id <- as.numeric(data$Subject)
data <- data[data$time<16*24,]
data <- lagfn(data, lagvars=c("time","conc"), id="Subject", time="time", lagfirst = NA)
head(data)
#>   Subject  Wt Apgar ApgarInd  time dose conc event id time.lag conc.lag
#> 1      42 2.8     9     >= 5  14.0   NA 13.3     1  1       NA       NA
#> 2      42 2.8     9     >= 5  95.5   NA 13.9     1  1     14.0     13.3
#> 3      28 3.2     9     >= 5   2.0   NA 16.9     1  2       NA       NA
#> 4      30 1.8     8     >= 5   6.3   NA 17.9     1  3       NA       NA
#> 5      30 1.8     8     >= 5 226.3   NA 16.5     1  3      6.3     17.9
#> 6      56 0.6     4      < 5  20.0   NA 18.8     1  4       NA       NA

mph <- coxph(Surv(time.lag,time,event)~I(conc.lag>0 & conc.lag<=20) + I(conc.lag>20 & conc.lag<=30)
 + I(conc.lag>30)+ cluster(id),,data=data)
summary(mph)
#> Call:
#> coxph(formula = Surv(time.lag, time, event) ~ I(conc.lag > 0 & 
#>     conc.lag <= 20) + I(conc.lag > 20 & conc.lag <= 30) + I(conc.lag > 
#>     30), data = data, cluster = id)
#> 
#>   n= 95, number of events= 95 
#>    (59 observations deleted due to missingness)
#> 
#>                                         coef exp(coef) se(coef) robust se     z
#> I(conc.lag > 0 & conc.lag <= 20)TRUE  0.6546    1.9243   0.3566    0.2914 2.246
#> I(conc.lag > 20 & conc.lag <= 30)TRUE 0.3741    1.4537   0.3142    0.2621 1.427
#> I(conc.lag > 30)TRUE                      NA        NA   0.0000    0.0000    NA
#>                                       Pr(>|z|)  
#> I(conc.lag > 0 & conc.lag <= 20)TRUE    0.0247 *
#> I(conc.lag > 20 & conc.lag <= 30)TRUE   0.1535  
#> I(conc.lag > 30)TRUE                        NA  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                       exp(coef) exp(-coef) lower .95 upper .95
#> I(conc.lag > 0 & conc.lag <= 20)TRUE      1.924     0.5197    1.0870     3.406
#> I(conc.lag > 20 & conc.lag <= 30)TRUE     1.454     0.6879    0.8697     2.430
#> I(conc.lag > 30)TRUE                         NA         NA        NA        NA
#> 
#> Concordance= 0.544  (se = 0.023 )
#> Likelihood ratio test= 3.53  on 2 df,   p=0.2
#> Wald test            = 5.19  on 2 df,   p=0.07
#> Score (logrank) test = 3.49  on 2 df,   p=0.2,   Robust = 6.28  p=0.04
#> 
#>   (Note: the likelihood ratio and score tests assume independence of
#>      observations within a cluster, the Wald and robust score tests do not).
data$weight <- iiw(mph,data,"id","time",TRUE)
head(data)
#>   Subject  Wt Apgar ApgarInd  time dose conc event id time.lag conc.lag
#> 1      42 2.8     9     >= 5  14.0   NA 13.3     1  1       NA       NA
#> 2      42 2.8     9     >= 5  95.5   NA 13.9     1  1     14.0     13.3
#> 3      28 3.2     9     >= 5   2.0   NA 16.9     1  2       NA       NA
#> 4      30 1.8     8     >= 5   6.3   NA 17.9     1  3       NA       NA
#> 5      30 1.8     8     >= 5 226.3   NA 16.5     1  3      6.3     17.9
#> 6      56 0.6     4      < 5  20.0   NA 18.8     1  4       NA       NA
#>      weight
#> 1 1.0000000
#> 2 0.5196749
#> 3 1.0000000
#> 4 1.0000000
#> 5 0.5196749
#> 6 1.0000000
```
