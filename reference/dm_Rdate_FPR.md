# Make FPR design matrix for dates with R format.

`dm_Rdate_FPR` creates design matrices for false positive rate
regressions; can also be used to standardize dates.

## Usage

``` r
dm_Rdate_FPR(Rdate, Y, effect = "fixed", num_knots_FPR = NULL)
```

## Arguments

- Rdate:

  a vector of dates of R format

- Y:

  binary case/control status; 1 for case; 0 for controls

- effect:

  The design matrix for "random" or "fixed" effect; Default is "fixed".
  When specified as "fixed", it produces standardized R-format dates
  using control's mean and standard deviation; When specified as
  "random", it produces `num_knots_FPR` columns of design matrix for
  thin-plate regression splines (TPRS) fitting. One needs both "fixed"
  and "random" in a FPR regression formula in `model_options` to enable
  TPRS fitting. For example, `model_options$likelihood$FPR_formula` can
  be  
    
  `~ AGECAT+HIV+dm_Rdate_FPR(ENRLDATE,Y,"fixed")+dm_Rdate_FPR(ENRLDATE,Y,"random",10)`  
    
  means FPR regression with intercept, main effects for 'AGECAT' and
  'HIV', and TPRS bases for 'ENRLDATE' using 10 knots placed at 10
  equal-probability-spaced sample quantiles.

- num_knots_FPR:

  number of knots for FPR regression; default is `NULL` to accommodate
  fixed effect specification.

## Value

Design matrix for FPR regression:

- `Z_FPR_ctrl` transformed design matrix for FPR regression for controls

- `Z_FPR_case` transformed design matrix for borrowing FPR regression
  from controls to cases. It is obtained using control-standardization,
  and square-root the following matrix (\\\Omega\\\]) with
  (\\j_1\\,\\j_2\\) element being
  \$\$\Omega\_{j_1j_2}=\\knots\_{j_1}-knots\_{j_2}\\^3\$\$.

## See also

[`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)
