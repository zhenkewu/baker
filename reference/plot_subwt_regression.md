# visualize the subclass weight regression with a continuous covariate

visualize the subclass weight regression with a continuous covariate

## Usage

``` r
plot_subwt_regression(
  DIR_NPLCM,
  stratum_bool,
  case = 0,
  slice = 1,
  truth = NULL,
  RES_NPLCM = NULL
)
```

## Arguments

- DIR_NPLCM:

  File path to the folder containing posterior samples

- stratum_bool:

  a vector of TRUE/FALSE with TRUE indicating the rows of subjects to
  include

- case:

  1 for plotting cases, 0 for plotting controls; default to 0.

- slice:

  integer; specifies which slice of bronze-standard data to visualize;
  Default to 1.

- truth:

  a list of truths computed from true parameters in simulations;
  elements: Eti, FPR, PR_case,TPR; All default to `NULL` in real data
  analyses. Currently only works for one slice of bronze-standard
  measurements (in a non-nested model).

  - truth_subwt matrix of \# of rows = \# of subjects, \# columns:
    number of true subclasses

- RES_NPLCM:

  pre-read res_nplcm; default to NULL.

## Value

A figure of subclass regression curves

## See also

Other visualization functions:
[`plot.nplcm()`](https://zhenkewu.com/baker/reference/plot.nplcm.md),
[`plot_BrS_panel()`](https://zhenkewu.com/baker/reference/plot_BrS_panel.md),
[`plot_SS_panel()`](https://zhenkewu.com/baker/reference/plot_SS_panel.md),
[`plot_check_common_pattern()`](https://zhenkewu.com/baker/reference/plot_check_common_pattern.md),
[`plot_check_pairwise_SLORD()`](https://zhenkewu.com/baker/reference/plot_check_pairwise_SLORD.md),
[`plot_etiology_regression()`](https://zhenkewu.com/baker/reference/plot_etiology_regression.md),
[`plot_etiology_strat()`](https://zhenkewu.com/baker/reference/plot_etiology_strat.md),
[`plot_panels()`](https://zhenkewu.com/baker/reference/plot_panels.md),
[`plot_pie_panel()`](https://zhenkewu.com/baker/reference/plot_pie_panel.md)
