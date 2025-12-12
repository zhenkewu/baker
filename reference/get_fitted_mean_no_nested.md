# get model fitted mean for conditional independence model

get model fitted mean for conditional independence model

## Usage

``` r
get_fitted_mean_no_nested(
  slice,
  res_nplcm,
  model_options,
  data_nplcm,
  clean_options
)
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

- clean_options:

  see
  [`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

## Value

a list with model fitted means
