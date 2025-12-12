# create regressor summation equation used in regression for FPR

`create_bugs_regressor_FPR` creates linear product of coefficients and a
row of design matrix used in regression

## Usage

``` r
create_bugs_regressor_FPR(
  n,
  dm_nm = "dm_FPR",
  b_nm = "b",
  ind_nm = "j",
  sub_ind_nm = "k"
)
```

## Arguments

- n:

  the length of coefficients

- dm_nm:

  name of design matrix; default `"dm_FPR"`

- b_nm:

  name of the coefficients; default `"b"`

- ind_nm:

  name of the coefficient iterator; default `"j"`

- sub_ind_nm:

  name of the subject iterator; default `"k"`

## Value

a character string with linear product form
