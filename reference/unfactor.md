# Convert factor to numeric without losing information on the label

Convert factor to numeric without losing information on the label

## Usage

``` r
unfactor(f)
```

## Arguments

- f:

  A factor

## Value

A numeric vector

## Examples

``` r
unfactor(factor(c("1","3","3"),levels=c("1","3")))
#> [1] 1 3 3
# contrast this to:
as.numeric(factor(c("1","3","3"),levels=c("1","3"))) 
#> [1] 1 2 2
```
