# Plot silver-standard (SS) panel

Plot silver-standard (SS) panel

## Usage

``` r
plot_SS_panel(
  slice,
  data_nplcm,
  model_options,
  clean_options,
  bugs.dat,
  res_nplcm,
  bg_color,
  select_latent = NULL,
  exact = TRUE,
  top_SS = 1,
  cexval = 1,
  srtval = 0,
  prior_shape = "interval"
)
```

## Arguments

- slice:

  the index of measurement slice for SS.

- data_nplcm:

  See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- model_options:

  See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- clean_options:

  See
  [`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

- bugs.dat:

  Data input for the model fitting.

- res_nplcm:

  See
  [`nplcm_read_folder()`](https://zhenkewu.com/baker/reference/nplcm_read_folder.md)

- bg_color:

  A list with names "BrS", "SS", "pie" to specify background colors

- select_latent:

  a vector of character strings representing latent status. It is used
  for just plotting a subset of latent status. For example, you can
  specify `select_latent = "HINF"` to plot all latent status information
  relevant to `"HINF"`.

- exact:

  Default is `TRUE` to use `select_latent` as exact names of causes. If
  you want to specify a name and plot all single or combo causes with
  that name, specify it to be `FALSE`.

- top_SS:

  Numerical value to specify the rightmost limit on the horizontal axis
  for the SS panel.

- cexval:

  Default is 1 - size of text of the SS percentages.

- srtval:

  Default is 0 - the direction of the text for the SS percentages.

- prior_shape:

  `interval` or `boxplot` - for how to represent prior/posteriors of the
  TPR/FPRs of measurements.

## Value

plotting function

## See also

Other visualization functions:
[`plot.nplcm()`](https://zhenkewu.com/baker/reference/plot.nplcm.md),
[`plot_BrS_panel()`](https://zhenkewu.com/baker/reference/plot_BrS_panel.md),
[`plot_check_common_pattern()`](https://zhenkewu.com/baker/reference/plot_check_common_pattern.md),
[`plot_check_pairwise_SLORD()`](https://zhenkewu.com/baker/reference/plot_check_pairwise_SLORD.md),
[`plot_etiology_regression()`](https://zhenkewu.com/baker/reference/plot_etiology_regression.md),
[`plot_etiology_strat()`](https://zhenkewu.com/baker/reference/plot_etiology_strat.md),
[`plot_panels()`](https://zhenkewu.com/baker/reference/plot_panels.md),
[`plot_pie_panel()`](https://zhenkewu.com/baker/reference/plot_pie_panel.md),
[`plot_subwt_regression()`](https://zhenkewu.com/baker/reference/plot_subwt_regression.md)
