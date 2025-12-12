# Set true positive rate (TPR) prior ranges for bronze-standard (BrS) data

`set_prior_tpr_BrS_NoNest` is for for conditional independence models.
We currently also use it for conditional dependence model: subclass TPRs
are independently assigned a beta prior. Ongoing work will enable
specifying priors for the marginal TPR, i.e. TPRs for a disease class,
not for the finer subclass.

## Usage

``` r
set_prior_tpr_BrS_NoNest(slice, model_options, data_nplcm)
```

## Arguments

- slice:

  The BrS measurement slice under consideration.

- model_options:

  See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)
  function.

- data_nplcm:

  See
  [`assign_model()`](https://zhenkewu.com/baker/reference/assign_model.md)
  function.

## Value

Parameters for the BrS dta TPR priors. It is a list of two lists (alpha
and beta). Alpha and beta are of the same length, the number of BrS
measurement slices. Each element of the alpha (beta) list is a numeric
vector for alpha (beta) parameters as in BETA distribution.

## See also

Other prior specification functions:
[`overall_uniform()`](https://zhenkewu.com/baker/reference/overall_uniform.md),
[`set_prior_tpr_SS()`](https://zhenkewu.com/baker/reference/set_prior_tpr_SS.md)
