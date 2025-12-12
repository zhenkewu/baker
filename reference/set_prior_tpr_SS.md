# Set true positive rate (TPR) prior ranges for silver-standard data.

Set true positive rate (TPR) prior ranges for silver-standard data.

## Usage

``` r
set_prior_tpr_SS(model_options, data_nplcm)
```

## Arguments

- model_options:

  See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)
  function.

- data_nplcm:

  See
  [`assign_model()`](https://zhenkewu.com/baker/reference/assign_model.md)
  function.

## Value

Parameters for the SS data TPR priors. It is a list of two lists (alpha
and beta). Alpha and beta are of the same length, the number of BrS
measurement slices. Each element of the alpha (beta) list is a numeric
vector for alpha (beta) parameters to specify Beta prior for TPRs.

## See also

Other prior specification functions:
[`overall_uniform()`](https://zhenkewu.com/baker/reference/overall_uniform.md),
[`set_prior_tpr_BrS_NoNest()`](https://zhenkewu.com/baker/reference/set_prior_tpr_BrS_NoNest.md)
