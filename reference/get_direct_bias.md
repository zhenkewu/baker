# Obtain direct bias that measure the discrepancy of a posterior distribution of pie and a true pie.

Obtain direct bias that measure the discrepancy of a posterior
distribution of pie and a true pie.

## Usage

``` r
get_direct_bias(DIR_list, truth = NULL, silent = FALSE)
```

## Arguments

- DIR_list:

  The list of where Bayesian results are stored

- truth:

  True etiologic fraction vector (must sum to 1) used to generate data;
  Default is `NULL`. If a vector is supplied, then only the first path
  in `DIR_LIST` is used.

- silent:

  Default is FALSE. To suppress printing messages, set to TRUE.

## Value

a list of length two. `diff` is the direct differences; `prb` is the
percent relative bias.
