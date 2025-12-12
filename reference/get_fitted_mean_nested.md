# get fitted mean for nested model with subclass mixing weights that are the same among cases

get fitted mean for nested model with subclass mixing weights that are
the same among cases

## Usage

``` r
get_fitted_mean_nested(
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

a matrix of no. of rows equal to retained MCMC samples, no. of columns
equal to the no. of measurement dimensions within a slice.
