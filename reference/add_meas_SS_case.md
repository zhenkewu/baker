# add likelihood for a SS measurement slice among cases (conditional independence)

add likelihood for a SS measurement slice among cases (conditional
independence)

## Usage

``` r
add_meas_SS_case(nslice, Mobs, prior, cause_list)
```

## Arguments

- nslice:

  the total number of SS measurement slices

- Mobs:

  see `data_nplcm` described in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- prior:

  see `model_options` described in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- cause_list:

  the list of causes in `model_options` described in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

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
[`add_meas_BrS_ctrl_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_param_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice.md),
[`add_meas_BrS_param_Nest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice_jags.md),
[`add_meas_BrS_param_Nest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice.md),
[`add_meas_BrS_param_NoNest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_subclass_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_subclass_Nest_Slice.md),
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
[`add_meas_BrS_ctrl_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_ctrl_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_param_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice.md),
[`add_meas_BrS_param_Nest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_Slice_jags.md),
[`add_meas_BrS_param_Nest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_Nest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice.md),
[`add_meas_BrS_param_NoNest_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_Slice_jags.md),
[`add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags()`](https://zhenkewu.com/baker/reference/add_meas_BrS_param_NoNest_reg_discrete_predictor_Slice_jags.md),
[`add_meas_BrS_subclass_Nest_Slice()`](https://zhenkewu.com/baker/reference/add_meas_BrS_subclass_Nest_Slice.md),
[`add_meas_SS_param()`](https://zhenkewu.com/baker/reference/add_meas_SS_param.md)
