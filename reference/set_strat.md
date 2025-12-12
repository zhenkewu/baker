# Stratification setup by covariates

`set_strat` makes group indicators based on `model_options$X_reg_*`

## Usage

``` r
set_strat(X, X_reg)
```

## Arguments

- X:

  A data frame of covariates

- X_reg:

  The vector of covariates that will stratify the analyses. These
  variables have to be categorical.

## Value

A list with following elements:

- `N_group` The number of groups

- `group` A vector of group indicator for every observation

## Details

the results from this function will help stratify etiology or FPR for
different strata; the ways of stratification for etiology and FPR can be
based on different covariates.
