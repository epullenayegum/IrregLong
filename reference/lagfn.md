# Create lagged versions the variables in data

Create lagged versions the variables in data

## Usage

``` r
lagfn(data, lagvars, id, time, lagfirst = NA)
```

## Arguments

- data:

  The data to be lagged

- lagvars:

  The names of the columns in the data to be lagged

- id:

  A character indicating which column of the data contains subject
  identifiers. ids are assumed to be consecutive integers, with the
  first subject having id 1

- time:

  A character indicating which column of the data contains the times at
  which each of the observations in data was made

- lagfirst:

  A vector giving the value of each lagged variable for the first time
  within each subject. This is helpful if, for example, time is the
  variable to be lagged and you know that all subjects entered the study
  at time zero

## Value

The original data frame with lagged variables added on as columns. For
example, if the data frame contains a variable named x giving the value
of x for each subject i at each visit j, the returned data frame will
contain a column named x.lag containing the value of x for subject i at
visit j-1. If j is the first visit for subject i, the lagged value is
set to NA

## Examples

``` r
library(nlme)
library(data.table)
data(Phenobarb)
head(Phenobarb)
#> Grouped Data: conc ~ time | Subject
#>   Subject  Wt Apgar ApgarInd time dose conc
#> 1       1 1.4     7     >= 5  0.0 25.0   NA
#> 2       1 1.4     7     >= 5  2.0   NA 17.3
#> 3       1 1.4     7     >= 5 12.5  3.5   NA
#> 4       1 1.4     7     >= 5 24.5  3.5   NA
#> 5       1 1.4     7     >= 5 37.0  3.5   NA
#> 6       1 1.4     7     >= 5 48.0  3.5   NA

data <- lagfn(Phenobarb,"time","Subject","time")
head(data)
#>   Subject  Wt Apgar ApgarInd time dose conc time.lag
#> 1      42 2.8     9     >= 5  0.0   28   NA       NA
#> 2      42 2.8     9     >= 5 12.0   28   NA      0.0
#> 3      42 2.8     9     >= 5 14.0   NA 13.3     12.0
#> 4      42 2.8     9     >= 5 23.7    7   NA     14.0
#> 5      42 2.8     9     >= 5 36.2    7   NA     23.7
#> 6      42 2.8     9     >= 5 47.8    7   NA     36.2
```
