# get marginal TPR and FPR for no nested model

get marginal TPR and FPR for no nested model

## Usage

``` r
get_marginal_rates_no_nested(slice, res_nplcm, model_options, data_nplcm)
```

## Arguments

- slice:

  the slice of BrS data that are modeled

- res_nplcm:

  matrix of MCMC samples

- model_options:

  see [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- data_nplcm:

  see [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

## Value

a matrix of no. of rows equal to retained MCMC samples, no. of columns
equal to the no. of measurement dimensions within a slice.
