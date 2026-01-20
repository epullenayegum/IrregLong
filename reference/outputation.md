# Create an outputted dataset for use with multiple outputation.

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
function creates a single outputted dataset.

## Usage

``` r
outputation(data, weights, singleobs, id, time, keep.first)
```

## Arguments

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

## Value

the outputted dataset.

## References

- Hoffman E, Sen P, Weinberg C. Within-cluster resampling. Biometrika
  2001; 88:1121-1134

- Follmann D, Proschan M, Leifer E. Multiple outputation: inference for
  complex clustered data by averaging analyses from independent data.
  Biometrics 2003; 59:420-429

- Pullenayegum EM. Multiple outputation for the analysis of longitudinal
  data subject to irregular observation. Statistics in Medicine (in
  press).

## See also

Other mo:
[`mo()`](https://epullenayegum.github.io/IrregLong/reference/mo.md)

## Examples

``` r
library(nlme)
data(Phenobarb)
library(survival)
library(data.table)
library(geepack)
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
head(data)
#> Grouped Data: conc ~ time | Subject
#>    Subject  Wt Apgar ApgarInd  time dose conc event id weight
#> 2        1 1.4     7     >= 5   2.0   NA 17.3     1 32  1.000
#> 12       1 1.4     7     >= 5 112.5   NA 31.0     1 32  0.478
#> 14       2 1.5     9     >= 5   2.0   NA  9.7     1 38  1.000
#> 20       2 1.5     9     >= 5  63.5   NA 24.6     1 38  1.000
#> 27       2 1.5     9     >= 5 135.5   NA 33.0     1 38  0.478
#> 29       3 1.5     6     >= 5   1.5   NA 18.0     1 11  1.000
data.output1 <-   outputation(data,data$weight,singleobs=FALSE,
id="id",time="time",keep.first=FALSE)
head(data.output1)
#> Grouped Data: conc ~ time | Subject
#>    Subject  Wt Apgar ApgarInd  time dose conc event id weight
#> 2        1 1.4     7     >= 5   2.0   NA 17.3     1 32  1.000
#> 12       1 1.4     7     >= 5 112.5   NA 31.0     1 32  0.478
#> 14       2 1.5     9     >= 5   2.0   NA  9.7     1 38  1.000
#> 20       2 1.5     9     >= 5  63.5   NA 24.6     1 38  1.000
#> 27       2 1.5     9     >= 5 135.5   NA 33.0     1 38  0.478
#> 29       3 1.5     6     >= 5   1.5   NA 18.0     1 11  1.000
data.output2 <-   outputation(data,data$weight,singleobs=FALSE,
id="id",time="time",keep.first=FALSE)
head(data.output2)
#> Grouped Data: conc ~ time | Subject
#>    Subject  Wt Apgar ApgarInd  time dose conc event id weight
#> 2        1 1.4     7     >= 5   2.0   NA 17.3     1 32  1.000
#> 14       2 1.5     9     >= 5   2.0   NA  9.7     1 38  1.000
#> 20       2 1.5     9     >= 5  63.5   NA 24.6     1 38  1.000
#> 27       2 1.5     9     >= 5 135.5   NA 33.0     1 38  0.478
#> 29       3 1.5     6     >= 5   1.5   NA 18.0     1 11  1.000
#> 36       3 1.5     6     >= 5  83.5   NA 23.8     1 11  1.000
data.output3 <-   outputation(data,data$weight,singleobs=FALSE,
id="id",time="time",keep.first=FALSE)
head(data.output3)
#> Grouped Data: conc ~ time | Subject
#>    Subject  Wt Apgar ApgarInd  time dose conc event id weight
#> 2        1 1.4     7     >= 5   2.0   NA 17.3     1 32  1.000
#> 12       1 1.4     7     >= 5 112.5   NA 31.0     1 32  0.478
#> 14       2 1.5     9     >= 5   2.0   NA  9.7     1 38  1.000
#> 20       2 1.5     9     >= 5  63.5   NA 24.6     1 38  1.000
#> 29       3 1.5     6     >= 5   1.5   NA 18.0     1 11  1.000
#> 36       3 1.5     6     >= 5  83.5   NA 23.8     1 11  1.000
# Note that the outputted dataset varies with each command run; outputation is done at random
```
