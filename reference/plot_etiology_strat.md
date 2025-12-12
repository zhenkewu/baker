# visualize the etiology estimates for each discrete levels

This function visualizes the etiology estimates against one discrete
covariate, e.g., age groups.

## Usage

``` r
plot_etiology_strat(
  DIR_NPLCM,
  strata_weights = "empirical",
  truth = NULL,
  RES_NPLCM = NULL,
  show_levels = 0,
  is_plot = TRUE,
  VERBOSE = TRUE
)
```

## Arguments

- DIR_NPLCM:

  File path to the folder containing posterior samples

- strata_weights:

  a vector of weights that sum to one; for each pathogen the weights
  specify how the j-th etiology fraction should be combined across all
  levels of the discrete predictors in the data; default is
  `"empirical"` to use empirical weights (observed fractions of subjects
  across strata).

- truth:

  a list of true values, e.g.,
  `truth=list(allEti = <a list of etiology fractions, each of identical length - the # of strata; >)`;
  if available, will be shown in thicker red solid vertical lines.

- RES_NPLCM:

  pre-read `res_nplcm`; default to `NULL`.

- show_levels:

  a vector of integers less than or equal to the total number of levels
  of strata; default to `0` for overall.

- is_plot:

  default to TRUE, plotting the figures; if `FALSE` only returning
  summaries

- VERBOSE:

  default to `TRUE`, print actual meanings of the levels

## Value

plotting function

## See also

Other visualization functions:
[`plot.nplcm()`](https://zhenkewu.com/baker/reference/plot.nplcm.md),
[`plot_BrS_panel()`](https://zhenkewu.com/baker/reference/plot_BrS_panel.md),
[`plot_SS_panel()`](https://zhenkewu.com/baker/reference/plot_SS_panel.md),
[`plot_check_common_pattern()`](https://zhenkewu.com/baker/reference/plot_check_common_pattern.md),
[`plot_check_pairwise_SLORD()`](https://zhenkewu.com/baker/reference/plot_check_pairwise_SLORD.md),
[`plot_etiology_regression()`](https://zhenkewu.com/baker/reference/plot_etiology_regression.md),
[`plot_panels()`](https://zhenkewu.com/baker/reference/plot_panels.md),
[`plot_pie_panel()`](https://zhenkewu.com/baker/reference/plot_pie_panel.md),
[`plot_subwt_regression()`](https://zhenkewu.com/baker/reference/plot_subwt_regression.md)
