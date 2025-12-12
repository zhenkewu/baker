# Obtain Integrated Squared Aitchison Distance, Squared Bias and Variance (both on Central Log-Ratio transformed scale) that measure the discrepancy of a posterior distribution of pie and a true pie.

The result is equivalent to Euclidean-type calculation after the
compositional vector (e.g., etiologic fraction) is centered-log-ratio
(CLRB) transformed. For simulation only.

## Usage

``` r
get_metric(DIR_NPLCM, truth)
```

## Arguments

- DIR_NPLCM:

  File path where Bayesian results are stored

- truth:

  True etiologic fraction vector (must sum to 1) used to generate data

## Value

a vector of (Integrated Squared Aitchison Distance (ISAD), bias-squared,
variance, truth)
