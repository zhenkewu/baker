# create regressor summation equation used in regression for etiology

`create_bugs_regressor_Eti` creates linear product of coefficients and a
row of design matrix used in regression

## Usage

``` r
create_bugs_regressor_Eti(
  n,
  dm_nm = "dm_Eti",
  b_nm = "betaEti",
  ind_nm = "j",
  sub_ind_nm = "k"
)
```

## Arguments

- n:

  the length of coefficients

- dm_nm:

  name of design matrix; default `"dm_Eti"`

- b_nm:

  name of the coefficients; default `"betaEti"`

- ind_nm:

  name of the coefficient iterator; default `"j"`

- sub_ind_nm:

  name of the subject iterator; default `"k"`

## Value

a character string with linear product form
