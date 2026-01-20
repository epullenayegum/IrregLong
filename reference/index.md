# Package index

## Examining the extent of irregularity

Functions for examining the extent of irregularity in observation times

- [`abacus.plot()`](https://epullenayegum.github.io/IrregLong/reference/abacus.plot.md)
  : Create an abacus plot Creates an abacus plot, depicting visits per
  subject over time
- [`extent.of.irregularity()`](https://epullenayegum.github.io/IrregLong/reference/extent.of.irregularity.md)
  : Measures of extent of visit irregularity Provides visual and numeric
  measures of the extent of irregularity in observation times in a
  longitudinal dataset

## Using inverse-intensity weighting

Using inverse-intensity weighting to account for irregular and
informative observation times

- [`iiw.weights()`](https://epullenayegum.github.io/IrregLong/reference/iiw.weights.md)
  : Compute inverse-intensity weights.
- [`iiwgee()`](https://epullenayegum.github.io/IrregLong/reference/iiwgee.md)
  : Fit an inverse-intensity weighted GEE.
- [`iiw()`](https://epullenayegum.github.io/IrregLong/reference/iiw.md)
  : Given a proportional hazards model for visit intensities, compute
  inverse-intensity weights.

## Multiple Outputation

Functions for performing multiple outputation

- [`mo()`](https://epullenayegum.github.io/IrregLong/reference/mo.md) :
  Multiple outputation for longitudinal data subject to irregular
  observation.
- [`outputation()`](https://epullenayegum.github.io/IrregLong/reference/outputation.md)
  : Create an outputted dataset for use with multiple outputation.

## Semi-parametric joint model

Functions for fitting the Liang semi-parametric joint model

- [`Liang()`](https://epullenayegum.github.io/IrregLong/reference/Liang.md)
  : Fit a semi-parametric joint model
- [`Liangint()`](https://epullenayegum.github.io/IrregLong/reference/Liangint.md)
  : Fit a semi-parametric joint model, incorporating intercept
  estimation
- [`create.bootstrapped.dataset()`](https://epullenayegum.github.io/IrregLong/reference/create.bootstrapped.dataset.md)
  : Create a single bootstrap sample for clustered data For clustered
  data, create a bootstrapped sample by sampling, with replacement, the
  same number of clusters as in the original dataset.

## Low-level functions

Functions for constructing the inverse-intensity weights manually

- [`addcensoredrows()`](https://epullenayegum.github.io/IrregLong/reference/addcensoredrows.md)
  : Add rows corresponding to censoring times to a longitudinal dataset
- [`lagfn()`](https://epullenayegum.github.io/IrregLong/reference/lagfn.md)
  : Create lagged versions the variables in data
