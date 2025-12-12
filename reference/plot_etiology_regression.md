# visualize the etiology regression with a continuous covariate

This function visualizes the etiology regression against one continuous
covariate, e.g., enrollment date. (NB: dealing with NoA,
multiple-pathogen causes, other continuous covariates? also there this
function only plots the first slice - so generalization may be useful -
give users an option to choose slice s; currently default to the first
slice.)

## Usage

``` r
plot_etiology_regression(
  DIR_NPLCM,
  stratum_bool,
  slice = 1,
  plot_basis = FALSE,
  truth = NULL,
  RES_NPLCM = NULL,
  do_plot = TRUE,
  do_rug = TRUE,
  return_metric = TRUE,
  plot_ma_dots = FALSE
)
```

## Arguments

- DIR_NPLCM:

  File path to the folder containing posterior samples

- stratum_bool:

  a vector of TRUE/FALSE with TRUE indicating the rows of subjects to
  include

- slice:

  integer; specifies which slice of bronze-standard data to visualize;
  Default to 1.

- plot_basis:

  TRUE for plotting basis functions; Default to FALSE

- truth:

  a list of truths computed from true parameters in simulations;
  elements: Eti, FPR, PR_case,TPR; All default to `NULL` in real data
  analyses. Currently only works for one slice of bronze-standard
  measurements (in a non-nested model).

  - Eti matrix of \# of rows = \# of subjects, \# columns:
    `length(cause_list)` for Eti

  - FPR matrix of \# of rows = \# of subjects, \# columns:
    `ncol(data_nplcm$Mobs$MBS$MBS1)`

  - PR_case matrix of \# of rows = \# of subjects, \# columns:
    `ncol(data_nplcm$Mobs$MBS$MBS1)`

  - TPR a vector of length identical to `PR_case`

- RES_NPLCM:

  pre-read res_nplcm; default to NULL.

- do_plot:

  TRUE for plotting

- do_rug:

  TRUE for plotting

- return_metric:

  TRUE for showing overall mean etiology, quantiles, s.d., and if
  `truth$Eti` is supplied, coverage, bias, truth and integrated mean
  squared errors (IMSE).

- plot_ma_dots:

  plot moving averages among case and controls if TRUE; Default to
  FALSE.

## Value

A figure of etiology regression curves and some marginal positive rate
assessment of model fit; See example for the legends.

## References

See example figures

- A Figure using simulated data for six pathogens:
  <https://github.com/zhenkewu/baker/blob/master/inst/figs/visualize_etiology_regression_SITE=1.pdf>

- The legends for the figure above:
  <https://github.com/zhenkewu/baker/blob/master/inst/figs/legends_visualize_etiology_regression.png>

## See also

Other visualization functions:
[`plot.nplcm()`](https://zhenkewu.com/baker/reference/plot.nplcm.md),
[`plot_BrS_panel()`](https://zhenkewu.com/baker/reference/plot_BrS_panel.md),
[`plot_SS_panel()`](https://zhenkewu.com/baker/reference/plot_SS_panel.md),
[`plot_check_common_pattern()`](https://zhenkewu.com/baker/reference/plot_check_common_pattern.md),
[`plot_check_pairwise_SLORD()`](https://zhenkewu.com/baker/reference/plot_check_pairwise_SLORD.md),
[`plot_etiology_strat()`](https://zhenkewu.com/baker/reference/plot_etiology_strat.md),
[`plot_panels()`](https://zhenkewu.com/baker/reference/plot_panels.md),
[`plot_pie_panel()`](https://zhenkewu.com/baker/reference/plot_pie_panel.md),
[`plot_subwt_regression()`](https://zhenkewu.com/baker/reference/plot_subwt_regression.md)
