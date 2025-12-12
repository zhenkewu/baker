# visualize trend of pathogen observation rate for NPPCR data (both cases and controls)

visualize trend of pathogen observation rate for NPPCR data (both cases
and controls)

## Usage

``` r
visualize_season(data_nplcm, patho, slice = 1, slice_SS = 1)
```

## Arguments

- data_nplcm:

  Data set produced by
  [`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

- patho:

  the index of pathogen

- slice:

  the slice of BrS data for visualization; default is 1.

- slice_SS:

  the slice of SS data to add onto BrS plots; default is 1, usually
  representing blood culture measurements.

## Value

A figure with smoothed positive rate and confidence bands for cases and
controls, respectively. The right margin shows marginal positive rates;
all rates are only among the subset of case subjects who had non-missing
responses for a measured agent (e.g., pathogen); similarly for controls.

## Details

This function shows observed positive rate for continuous
covariates,e.g., enrollment date in PERCH application. Smoothing is done
by penalized splines implemented by `mgcv` package. The penalized spline
smoothing term is constructed by
[`mgcv::smooth.construct.ps.smooth.spec()`](https://rdrr.io/pkg/mgcv/man/smooth.construct.ps.smooth.spec.html).
The window size of the moving averages currently is set to 60 days.

## See also

Other exploratory data analysis functions:
[`get_top_pattern()`](https://zhenkewu.com/baker/reference/get_top_pattern.md),
[`plot_logORmat()`](https://zhenkewu.com/baker/reference/plot_logORmat.md),
[`show_individual()`](https://zhenkewu.com/baker/reference/show_individual.md),
[`summarize_BrS()`](https://zhenkewu.com/baker/reference/summarize_BrS.md),
[`summarize_SS()`](https://zhenkewu.com/baker/reference/summarize_SS.md)
