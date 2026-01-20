# Add rows corresponding to censoring times to a longitudinal dataset

Add rows corresponding to censoring times to a longitudinal dataset

## Usage

``` r
addcensoredrows(data, maxfu, tinvarcols, id, time, event)
```

## Arguments

- data:

  The dataset to which rows are to be added. The data should have one
  row per observation

- maxfu:

  The maximum follow-up time per subject. If all subjects have the same
  follow-up time, this can be supplied as a single number. Otherwise,
  maxfu should be a dataframe with the first column specifying subject
  identifiers and the second giving the follow-up time for each subject.

- tinvarcols:

  A vector of column numbers corresponding to variables in data that are
  time-invariant.

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

## Value

The original dataset with extra rows corresponding to censoring times

## Examples

``` r
x <- c(1:3,1:2,1:5)
x0 <- c(rep(2,3),rep(0,2),rep(1,5))
id <- c(rep(1,3),rep(2,2),rep(3,5))
time <- c(0,4,6,2,3,1,3,5,6,7)
event <- c(1,1,1,0,1,0,1,1,1,1)
data <- as.data.frame(cbind(x,id,time,event,x0))
addcensoredrows(data,maxfu=8,id="id",time="time",tinvarcols=5,event="event")
#>     x id time event x0
#> 1   1  1    0     1  2
#> 2   2  1    4     1  2
#> 3   3  1    6     1  2
#> 11 NA  1    8     0  2
#> 4   1  2    2     0  0
#> 5   2  2    3     1  0
#> 12 NA  2    8     0  0
#> 6   1  3    1     0  1
#> 7   2  3    3     1  1
#> 8   3  3    5     1  1
#> 9   4  3    6     1  1
#> 10  5  3    7     1  1
#> 13 NA  3    8     0  1


x <- c(1:3,1:2,1:5)
x0 <- c(rep(2,3),rep(0,2),rep(1,5))
id <- c(rep(1,3),rep(2,2),rep(3,5))
time <- c(0,4,6,2,3,1,3,5,6,7)
event <- c(1,1,1,0,1,0,1,1,1,1)
data <- as.data.frame(cbind(x,id,time,event,x0))
maxfu.id <- 1:3
maxfu.time <- c(6,5,8)
maxfu <- cbind(maxfu.id,maxfu.time)
maxfu <- as.data.frame(maxfu)
addcensoredrows(data,maxfu=maxfu,id="id",time="time",tinvarcols=5,event="event")
#>     x id time event x0
#> 1   1  1    0     1  2
#> 2   2  1    4     1  2
#> 3   3  1    6     1  2
#> 4   1  2    2     0  0
#> 5   2  2    3     1  0
#> 11 NA  2    5     0  0
#> 6   1  3    1     0  1
#> 7   2  3    3     1  1
#> 8   3  3    5     1  1
#> 9   4  3    6     1  1
#> 10  5  3    7     1  1
#> 12 NA  3    8     0  1
```
