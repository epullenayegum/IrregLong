# Multiple outputation for longitudinal data subject to irregular observation.

Multiple outputation is a procedure whereby excess observations are
repeatedly randomly sampled and discarded. The method was originally
developed to handle clustered data where cluster size is informative,
for example when studying pups in a litter. In this case, analysis that
ignores cluster size results in larger litters being over-represented in
a marginal analysis. Multiple outputation circumvents this problem by
randomly selecting one observation per cluster. Multiple outputation has
been further adapted to handle longitudinal data subject to irregular
observation; here the probability of being retained on any given
outputation is inversely proportional to the visit intensity. This
function creates multiply outputted datasets, analyses each separately,
and combines the results to produce a single estimate.

## Usage

``` r
mo(
  noutput,
  fn,
  data,
  weights,
  singleobs,
  id,
  time,
  keep.first,
  var = TRUE,
  ...
)
```

## Arguments

- noutput:

  the number of outputations to be used

- fn:

  the function to be applied to the outputted datasets. fn should return
  a vector or scalar; if var=TRUE the second column of fn should be an
  estimate of standard error.

- data:

  the original dataset on which multiple outputation is to be performed

- weights:

  the weights to be used in the outputation, i.e. the inverse of the
  probability that a given observation will be selected in creating an
  outputted dataset. Ignored if singleobs=TRUE

- singleobs:

  logical variable indicating whether a single observation should be
  retained for each subject

- id:

  character string indicating which column of the data identifies
  subjects

- time:

  character string indicating which column of the data contains the time
  at which the visit occurred

- keep.first:

  logical variable indicating whether the first observation should be
  retained with probability 1. This is useful if the data consists of an
  observation at baseline followed by follow-up at stochastic time
  points.

- var:

  logical variable indicating whether fn returns variances in addition
  to point estimates

- ...:

  other arguments to fn.

## Value

a list containing the multiple outputation estimate of the function fn
applied to the data, its standard error, and the relative efficiency of
using noutput outputations as opposed to an infinite number

## References

- Hoffman E, Sen P, Weinberg C. Within-cluster resampling. Biometrika
  2001; 88:1121-1134

- Follmann D, Proschan M, Leifer E. Multiple outputation: inference for
  complex clustered data by averaging analyses from independent data.
  Biometrics 2003; 59:420-429

- Pullenayegum EM. Multiple outputation for the analysis of longitudinal
  data subject to irregular observation. Statistics in Medicine (in
  press)

.

## See also

Other mo:
[`outputation()`](https://epullenayegum.github.io/IrregLong/reference/outputation.md)

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
id="id",time="time",event="event",data=data, invariant="id",lagvars=c("time","conc"),maxfu=16*24,
         lagfirst=c(0,0),first=TRUE)
wt <- i$iiw.weight
wt[wt>quantile(i$iiw.weight,0.95)] <- quantile(i$iiw.weight,0.95)
data$wt <- wt
reg <- function(data){
est <- summary(geeglm(conc~I(time^3) + log(time), id=id,data=data))$coefficients[,1:2]
est <- data.matrix(est)
return(est)
}

mo(20,reg,data,wt,singleobs=FALSE,id="id",time="time",keep.first=FALSE)
#> $est
#> [1]  1.68e+01 -6.60e-07  3.11e+00
#> 
#> $se
#> [1] 1.103673 0.000343 0.670909
#> 
#> $RE.MO
#> [1] 1.01 1.00 1.01
#> 
# On outputation, the dataset contains small numbers of observations per subject
# and hence the GEE sandwich variance estimate underestimates variance; this is why
# the outputation-based variance estimate fails. This can be remedied by using a
# sandwich variance error correction (e.g. Fay-Graubard, Mancl-DeRouen).
```
