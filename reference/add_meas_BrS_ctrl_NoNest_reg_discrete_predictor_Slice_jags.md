# add a likelihood component for a BrS measurement slice among controls

regression model without nested subclasses; discrete

## Usage

``` r
add_meas_BrS_ctrl_NoNest_reg_discrete_predictor_Slice_jags(
  s,
  Mobs,
  cause_list,
  ppd = NULL
)
```

## Arguments

- s:

  the slice

- Mobs:

  See `data_nplcm` described in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- cause_list:

  the list of causes in `data_nplcm` described in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- ppd:

  Default is NULL; Set to TRUE for enabling posterior predictive
  checking.

## Value

a list of two elements: the first is `plug`, the .bug code; the second
is `parameters` that stores model parameters introduced by this plugged
measurement slice

## See also

Other likelihood specification functions:
[`add_meas_BrS_case_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_Nest_Slice.md),
[`add_meas_BrS_case_Nest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_Nest_Slice_jags.md),
[`add_meas_BrS_case_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_Slice.md),
[`add_meas_BrS_case_NoNest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_Slice_jags.md),
[`add_meas_BrS_case_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_case_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_ctrl_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_Nest_Slice.md),
[`add_meas_BrS_ctrl_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_NoNest_Slice.md),
[`add_meas_BrS_ctrl_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_param_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice.md),
[`add_meas_BrS_param_Nest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice_jags.md),
[`add_meas_BrS_param_Nest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice.md),
[`add_meas_BrS_param_NoNest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_subclass_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_subclass_Nest_Slice.md),
[`add_meas_SS_case()`](https://zhenkewu.com/baker/reference/add_meas_SS_case.md),
[`add_meas_SS_param()`](https://zhenkewu.com/baker/reference/add_meas_SS_param.md)

Other plug-and-play functions:
[`add_meas_BrS_case_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_Nest_Slice.md),
[`add_meas_BrS_case_Nest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_Nest_Slice_jags.md),
[`add_meas_BrS_case_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_Slice.md),
[`add_meas_BrS_case_NoNest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_Slice_jags.md),
[`add_meas_BrS_case_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_case_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_case_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_ctrl_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_Nest_Slice.md),
[`add_meas_BrS_ctrl_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_NoNest_Slice.md),
[`add_meas_BrS_ctrl_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_param_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice.md),
[`add_meas_BrS_param_Nest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice_jags.md),
[`add_meas_BrS_param_Nest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice.md),
[`add_meas_BrS_param_NoNest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_subclass_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_subclass_Nest_Slice.md),
[`add_meas_SS_case()`](https://zhenkewu.com/baker/reference/add_meas_SS_case.md),
[`add_meas_SS_param()`](https://zhenkewu.com/baker/reference/add_meas_SS_param.md)
