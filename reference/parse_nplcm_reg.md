# parse regression components (either false positive rate or etiology regression) for fitting npLCM; Only use this when formula is not `NULL`.

parse regression components (either false positive rate or etiology
regression) for fitting npLCM; Only use this when formula is not `NULL`.

## Usage

``` r
parse_nplcm_reg(form, data_nplcm, silent = TRUE)
```

## Arguments

- form:

  regression formula

- data_nplcm:

  data object for
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md); may
  contain covariates X; must have case-control status Y.

- silent:

  Default is `TRUE` for no message about covariates; `FALSE` otherwise.

## Value

`TRUE` for doing regression; `FALSE` otherwise.
