# Plot three-panel figures for nested partially-latent model results

`plot_panels()` visualizes the model outputs for communicating how the
data inform final latent disease status (etiology). It works for
singleton or combo etiologies.

## Usage

``` r
plot_panels(
  DIR_NPLCM,
  slices = "all",
  bg_color = list(BrS = "lavenderblush", SS = "mistyrose", pie = "antiquewhite"),
  select_latent = NULL,
  exact = TRUE,
  SS_upperlimit = 1,
  eti_upperlimit = 1,
  silent = TRUE,
  ref_eti0 = NULL,
  is_plot = TRUE
)
```

## Arguments

- DIR_NPLCM:

  File path to the folder containing posterior samples

- slices:

  DEFAULT is "all" - to plot all measurements; Otherwise, one can
  specify a list: `list(MBS=c(1,3),MSS=1)` means to plot the 1st and 3rd
  slice of BrS measurements and 1st of SS measurement.

- bg_color:

  A list with names "BrS", "SS", "pie" to specify background colors. The
  current default is
  `list(BrS = "lavenderblush", SS = "mistyrose", pie="antiquewhite")`.
  If no background is intended, specify as NULL or for a particular
  measurement, e.g., `BrS = NULL`.

- select_latent:

  a vector of character strings representing latent status. It is used
  for just plotting a subset of latent status. For example, you can
  specify `select_latent = "HINF"` to plot all latent status information
  relevant to `"HINF"`.

- exact:

  Default is `TRUE` to use `select_latent` as exact names of causes. If
  you want to specify a name and plot all single or combo causes with
  that name, specify it to be `FALSE`.

- SS_upperlimit:

  The upper limit of horizontal bar for the silver-standard subpanel
  (the middle panel). The default value is .25.

- eti_upperlimit:

  The upper limit of horizontal bar for the etiology posterior subpanel
  (the rightmost panel). The default value is .4

- silent:

  Default is TRUE to not print any warning messages; FALSE otherwise.

- ref_eti0:

  reference quantiles and means; a list: pEti_ref_q, pEti_ref_mean_ord

- is_plot:

  default to `TRUE` for plotting only; set to `FALSE` if to get summary.

## Value

A figure with two or three columns (if `is_plot=TRUE`); otherwise, it
provide posterior summaries of Etiology information to used by
[`print.summary.nplcm.no_reg()`](https://zhenkewu.com/baker/reference/print.summary.nplcm.no_reg.md)

## Details

Missing data for BrS or SS are dropped when calculating observed
measurement positive rates

## See also

Other visualization functions:
[`plot.nplcm()`](https://zhenkewu.com/baker/reference/plot.nplcm.md),
[`plot_BrS_panel()`](https://zhenkewu.com/baker/reference/plot_BrS_panel.md),
[`plot_SS_panel()`](https://zhenkewu.com/baker/reference/plot_SS_panel.md),
[`plot_check_common_pattern()`](https://zhenkewu.com/baker/reference/plot_check_common_pattern.md),
[`plot_check_pairwise_SLORD()`](https://zhenkewu.com/baker/reference/plot_check_pairwise_SLORD.md),
[`plot_etiology_regression()`](https://zhenkewu.com/baker/reference/plot_etiology_regression.md),
[`plot_etiology_strat()`](https://zhenkewu.com/baker/reference/plot_etiology_strat.md),
[`plot_pie_panel()`](https://zhenkewu.com/baker/reference/plot_pie_panel.md),
[`plot_subwt_regression()`](https://zhenkewu.com/baker/reference/plot_subwt_regression.md)
