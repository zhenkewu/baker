# softmax

uses logsumexp trick to prevent numerical overflow

## Usage

``` r
softmax(x)
```

## Arguments

- x:

  a vector of numbers

## Value

a vector of positive values that sum to one.

## Examples

``` r
softmax2 <- function(x) exp(x) / sum(exp(x))
softmax(c(1, 2, 3) * 1000)  # NaN NaN NaN
#> [1] 0 0 1
softmax2(c(1, 2, 3) * 1000)  # 0 0 1
#> [1] NaN NaN NaN
```
