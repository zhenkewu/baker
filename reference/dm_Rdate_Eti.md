# Make etiology design matrix for dates with R format.

`dm_Rdate_Eti` creates design matrices for etiology regressions.

## Usage

``` r
dm_Rdate_Eti(Rdate, Y, num_knots_Eti, basis_Eti = "ncs")
```

## Arguments

- Rdate:

  a vector of dates of R format

- Y:

  binary case/control status; 1 for case; 0 for controls

- num_knots_Eti:

  number of knots for etiology regression

- basis_Eti:

  the type of basis functions to use for etiology regression. It can be
  "ncs" (natural cubic splines) or "tprs" (thin-plate regression
  splines). Default is "ncs". "tprs" will be implemented later.

## Value

Design matrix for etiology regression:

- `Z_Eti` transformed design matrix for etiology regression

## Details

It is used in `model_options$likeihood$Eti_formula`. For example, one
can specify it as:  
  
`~ AGECAT+HIV+dm_Rdate_Eti(ENRLDATE,Y,5)`  
  
to call an etiology regression with intercept, main effects for 'AGECAT'
and 'HIV', and natural cubic spline bases for 'ENRLDATE' using 5 knots
defined as 5 equal-probability-spaced sample quantiles.

## See also

[`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)
