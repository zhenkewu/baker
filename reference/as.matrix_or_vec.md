# convert one column data frame to a vector

convert one column data frame to a vector

## Usage

``` r
as.matrix_or_vec(x)
```

## Arguments

- x:

  an one-column data.frame

## Value

a vector

## Details

JAGS cannot accept a data frame with one column; This function converts
it to a vector, which JAGS will allow.
