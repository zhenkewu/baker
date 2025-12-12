# visualize the PERCH etiology regression with a continuous covariate

This function is specifically designed for PERCH data, e.g., (NB:
dealing with NoA, multiple-pathogen causes, other continuous covariates?
also there this function only plots the first slice - so generalization
may be useful - give users an option to choose slice s; currently
default to the first slice.)

## Usage

``` r
plot_case_study(
  DIR_NPLCM,
  stratum_bool = stratum_bool,
  bugs.dat = NULL,
  slice = 1,
  RES_NPLCM = NULL,
  do_plot = TRUE,
  do_rug = FALSE,
  return_metric = TRUE
)
```

## Arguments

- DIR_NPLCM:

  File path to the folder containing posterior samples

- stratum_bool:

  integer; for this function, indicates which strata to plot

- bugs.dat:

  The posterior samples (loaded into the environment to save time) -\>
  default is NULL

- slice:

  integer; specifies which slice of bronze-standard data to visualize;
  Default to 1.

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

## Value

A figure of etiology regression curves and some marginal positive rate
assessment of model fit; See example for the legends.
