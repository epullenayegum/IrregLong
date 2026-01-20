# Compute inverse-intensity weights.

Since the vector of weights is ordered on id and time, if you intend to
merge these weights onto your original dataset it is highly recommended
that you sort the data before running iiw.weights. The loess.stabilize
option is designed to avoid time trends in the weights. This option fits
a loess smooth of the weights vs. time, then divides by the predicted
value.

## Usage

``` r
iiw.weights(
  formulaph,
  formulanull = NULL,
  data,
  id,
  time,
  event,
  lagvars,
  invariant = NULL,
  maxfu,
  lagfirst = lagfirst,
  first,
  stabilize.loess = FALSE
)
```

## Arguments

- formulaph:

  the formula for the proportional hazards model for the visit intensity
  that will be used to derive inverse-intensity weights. The formula
  should usually use the counting process format (i.e.
  Surv(start,stop,event)). If a frailty model is used, the cluster(id)
  term should appear before other covariates

- formulanull:

  if stabilised weights are to be used, the formula for the null model
  used to stabilise the weights

- data:

  data frame containing the variables in the model

- id:

  character string indicating which column of the data identifies
  subjects

- time:

  character string indicating which column of the data contains the time
  at which the visit occurred

- event:

  character string indicating which column of the data indicates whether
  or not a visit occurred. If every row corresponds to a visit, then
  this column will consist entirely of ones

- lagvars:

  a vector of variable names corresponding to variables which need to be
  lagged by one visit to fit the visit intensity model. Typically time
  will be one of these variables. The function will internally add
  columns to the data containing the values of the lagged variables from
  the previous visit. Values of lagged variables for a subject's first
  visit will be set to NA. To access these variables in specifying the
  proportional hazards formulae, add ".lag" to the variable you wish to
  lag. For example, if time is the variable for time, time.lag is the
  time of the previous visit

- invariant:

  a vector of variable names corresponding to variables in data that are
  time-invariant. It is not necessary to list every such variable, just
  those that are invariant and also included in the proportional hazards
  model

- maxfu:

  the maximum follow-up time(s). If everyone is followed for the same
  length of time, this can be given as a single value. Otherwise, maxfu
  should be a dataframe with the first column specifying subject
  identifiers and the second giving the follow-up time for each subject.

- lagfirst:

  A vector giving the value of each lagged variable for the first time
  within each subject. This is helpful if, for example, time is the
  variable to be lagged and you know that all subjects entered the study
  at time zero

- first:

  logical variable. If TRUE, the first observation for each individual
  is assigned an intensity of 1. This is appropriate if the first visit
  is a baseline visit at which recruitment to the study occurred; in
  this case the baseline visit is observed with probability 1.

- stabilize.loess:

  logical variable. If TRUE, additional stabilization is done by fitting
  a loess of the (stabilized) weights versus time, then dividing the
  observed weights by the predicted values

## Value

a vector of inverse-intensity weights, ordered on id then time

## Details

Given longitudinal data with irregular visit times, fit a Cox
proportional hazards model for the visit intensity, then use it to
compute inverse-intensity weights

## References

- Lin H, Scharfstein DO, Rosenheck RA. Analysis of Longitudinal data
  with Irregular, Informative Follow-up. Journal of the Royal
  Statistical Society, Series B (2004), 66:791-813

- Buzkova P, Lumley T. Longitudinal data analysis for generalized linear
  models with follow-up dependent on outcome-related variables. The
  Canadian Journal of Statistics 2007; 35:485-500.

## See also

Other iiw:
[`iiw()`](https://epullenayegum.github.io/IrregLong/reference/iiw.md),
[`iiwgee()`](https://epullenayegum.github.io/IrregLong/reference/iiwgee.md)

Other iiw:
[`iiw()`](https://epullenayegum.github.io/IrregLong/reference/iiw.md),
[`iiwgee()`](https://epullenayegum.github.io/IrregLong/reference/iiwgee.md)

## Examples

``` r
library(nlme)
data(Phenobarb)
library(survival)
library(geepack)
library(data.table)
Phenobarb$event <- 1-as.numeric(is.na(Phenobarb$conc))
data <- Phenobarb
data <- data[data$event==1,]
data$id <- as.numeric(data$Subject)
data <- data[data$time<16*24,]
i <- iiw.weights(Surv(time.lag,time,event)~I(conc.lag>0 & conc.lag<=20) +
I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30)+ cluster(id),
id="id",time="time",event="event",data=data,
invariant="id",lagvars=c("time","conc"),maxfu=16*24,lagfirst=0,first=TRUE)
data$weight <- i$iiw.weight
summary(i$m)
#> Call:
#> coxph(formula = Surv(time.lag, time, event) ~ I(conc.lag > 0 & 
#>     conc.lag <= 20) + I(conc.lag > 20 & conc.lag <= 30) + I(conc.lag > 
#>     30), data = datacox, cluster = id)
#> 
#>   n= 154, number of events= 95 
#>    (59 observations deleted due to missingness)
#> 
#>                                         coef exp(coef) se(coef) robust se     z
#> I(conc.lag > 0 & conc.lag <= 20)TRUE  0.7379    2.0915   0.3398    0.3066 2.406
#> I(conc.lag > 20 & conc.lag <= 30)TRUE 0.3261    1.3856   0.3004    0.2851 1.144
#> I(conc.lag > 30)TRUE                      NA        NA   0.0000    0.0000    NA
#>                                       Pr(>|z|)  
#> I(conc.lag > 0 & conc.lag <= 20)TRUE    0.0161 *
#> I(conc.lag > 20 & conc.lag <= 30)TRUE   0.2527  
#> I(conc.lag > 30)TRUE                        NA  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#>                                       exp(coef) exp(-coef) lower .95 upper .95
#> I(conc.lag > 0 & conc.lag <= 20)TRUE      2.091     0.4781    1.1467     3.815
#> I(conc.lag > 20 & conc.lag <= 30)TRUE     1.386     0.7217    0.7924     2.423
#> I(conc.lag > 30)TRUE                         NA         NA        NA        NA
#> 
#> Concordance= 0.559  (se = 0.022 )
#> Likelihood ratio test= 5.12  on 2 df,   p=0.08
#> Wald test            = 6.69  on 2 df,   p=0.04
#> Score (logrank) test = 5.31  on 2 df,   p=0.07,   Robust = 7.18  p=0.03
#> 
#>   (Note: the likelihood ratio and score tests assume independence of
#>      observations within a cluster, the Wald and robust score tests do not).
# can use to fit a weighted GEE
mw <- geeglm(conc ~ I(time^3) + log(time) , id=Subject, data=data, weights=weight)
summary(mw)
#> 
#> Call:
#> geeglm(formula = conc ~ I(time^3) + log(time), data = data, weights = weight, 
#>     id = Subject)
#> 
#>  Coefficients:
#>               Estimate    Std.err   Wald Pr(>|W|)    
#> (Intercept)  1.684e+01  1.049e+00 257.82   <2e-16 ***
#> I(time^3)   -6.410e-07  6.112e-08 109.99   <2e-16 ***
#> log(time)    3.064e+00  3.588e-01  72.89   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Correlation structure = independence 
#> Estimated Scale Parameters:
#> 
#>             Estimate Std.err
#> (Intercept)    47.31   14.66
#> Number of clusters:   59  Maximum cluster size: 6 
# agrees with results through the single command iiwgee
miiwgee <- iiwgee(conc ~ I(time^3) + log(time),
Surv(time.lag,time,event)~I(conc.lag>0 & conc.lag<=20) +
I(conc.lag>20 & conc.lag<=30) + I(conc.lag>30)+ cluster(id),
id="id",time="time",event="event",data=data,
invariant="id",lagvars=c("time","conc"),maxfu=16*24,lagfirst=0,first=TRUE)
summary(miiwgee$geefit)
#> 
#> Call:
#> geeglm(formula = formulagee, family = family, data = data, weights = useweight, 
#>     id = iddup, corstr = "independence")
#> 
#>  Coefficients:
#>              Estimate   Std.err  Wald Pr(>|W|)    
#> (Intercept)  1.65e+01  1.02e+00 262.0   <2e-16 ***
#> I(time^3)   -6.63e-07  7.65e-08  75.2   <2e-16 ***
#> log(time)    3.25e+00  4.05e-01  64.4    1e-15 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Correlation structure = independence 
#> Estimated Scale Parameters:
#> 
#>             Estimate Std.err
#> (Intercept)     49.1      17
#> Number of clusters:   59  Maximum cluster size: 6 
```
