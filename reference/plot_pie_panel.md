# Plot etiology (pie) panel

Plot etiology (pie) panel

## Usage

``` r
plot_pie_panel(
  model_options,
  res_nplcm,
  bugs.dat,
  bg_color,
  select_latent = NULL,
  exact = TRUE,
  top_pie = 1,
  label_size = 1,
  ref_eti = NULL,
  is_plot = TRUE
)
```

## Arguments

- model_options:

  See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- res_nplcm:

  See
  [`nplcm_read_folder()`](https://zhenkewu.com/baker/reference/nplcm_read_folder.md)

- bugs.dat:

  Data input for the model fitting.

- bg_color:

  A list with names "BrS", "SS", "pie" to specify background colors

- select_latent:

  a vector of character strings representing latent status. It is used
  for just plotting a subset of latent status. For example, you can
  specify `select_latent = "HINF"`

- exact:

  Default is `TRUE` to use `select_latent` as exact names of causes. If
  you want to specify a name and plot all single or combo causes with
  that name, specify it to be `FALSE`. to plot all latent status
  information relevant to `"HINF"`.

- top_pie:

  Numerical value to specify the rightmost limit on the horizontal axis
  for the pie panel.

- label_size:

  the size of latent status labels on the right margin

- ref_eti:

  reference quantiles and means; a list: pEti_ref_q, pEti_ref_mean_ord

- is_plot:

  default to `TRUE` for plotting only; set to `FALSE` if to get summary.

## Value

plotting function.

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
[`plot_subwt_regression()`](https://zhenkewu.com/baker/reference/plot_subwt_regression.md)
