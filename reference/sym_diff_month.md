# get symmetric difference of months from two vector of R-format dates

`sym_diff_month` evaluates the symmetric difference between two sets of
R-formatted date

## Usage

``` r
sym_diff_month(Rdate1, Rdate2)
```

## Arguments

- Rdate1, Rdate2:

  R-formatted R dates. See
  [`as.Date()`](https://rdrr.io/r/base/as.Date.html)

## Value

`NULL` if no difference; the set of different months otherwise.
