# Changelog

## IrregLong 0.4.1

## IrregLong 0.3.4

CRAN release: 2022-02-21

### Minor changes

- Added error message to extent.of.irregularity for when formula and/or
  ncutpts are not specified in the absence of scheduled visit times
- made use of data.table package when lagging variables and adding
  censored rows

### Bug fixes

- extent.of.irregularity previously only allowed the time variable to be
  called “time”; this has now been fixed

## IrregLong 0.3.3

CRAN release: 2020-08-24

### Minor changes

- due to archival of the frailtypack package, frailty models are now
  estimated using coxph rather than frailtyPenal

## IrregLong 0.3.2

CRAN release: 2020-06-24

### Bug fixes

- iiw.weights failed to consider the interval between last observed
  visit and end of follow-up time when maxfu was supplied as a single
  number. This has been fixed.

## IrregLong 0.3.1

CRAN release: 2020-05-23

### Bug fixes

- Fixed a coding error that led extent.of.irregularity to fail when
  scheduledtimes was provided.

## IrregLong 0.3.0

CRAN release: 2020-05-01

### Major changes

- Added function (extent.of.irregularity) to provide graphical and
  numeric summaries of the extent of irregularity.

### Minor changes

- Added an option in the Liang function to provide as an argument a
  function computing the value of the covariates X for any subject at
  any given time t, not just the observation times.

## IrregLong 0.2.0

CRAN release: 2020-01-24

### Major changes

- When using lagged variable, there is an option to include a value for
  the first lagged observation within each subject.
- A function (abacus.plot) for producing an abacus plot has been added.
- The Liang model now allows for time-invariant covariates to be
  included in the visit process model.

### Bug fixes

- The documentation (examples and vignette) has been updated to use
  dosing by weight in the measurement frequency model.
- tinvarcols can be omitted from addcensoredrows, iiwgee, and
  iiw.weights if maxfu=NULL.
- lagfn and addcensoredrows no longer require subject identifiers to be
  numeric.
