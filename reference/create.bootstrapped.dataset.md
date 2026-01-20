# Create a single bootstrap sample for clustered data For clustered data, create a bootstrapped sample by sampling, with replacement, the same number of clusters as in the original dataset.

Create a single bootstrap sample for clustered data For clustered data,
create a bootstrapped sample by sampling, with replacement, the same
number of clusters as in the original dataset.

## Usage

``` r
create.bootstrapped.dataset(data, idname)
```

## Arguments

- data:

  The dataset to be resampled

- idname:

  A character indicating which column of the data contains subject
  identifiers.

## Details

This function is designed to assist in computing bootstrap standard
errors when working with longitudinal data. Given longitudinal data with
multiple rows per subject, it will sample subjects, with replacement, n
times, where n is the number of subjects in the original dataset. In the
bootstrapped dataset, each resample has its own id.
