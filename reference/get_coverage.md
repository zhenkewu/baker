# Obtain coverage status from a result folder

Obtain coverage status from a result folder

## Usage

``` r
get_coverage(DIR_NPLCM, truth)
```

## Arguments

- DIR_NPLCM:

  Path to where Bayesian results are stored

- truth:

  True etiologic fraction vector (must sum to 1) used to generate data.

## Value

A logic vector of length as `truth`. 1 for covered; 0 for not.
