# Pick parameters in the Beta distribution to match the specified range

`beta_parms_from_quantiles` produces prior Beta parameters for the true
positive rates (TPR)

## Usage

``` r
beta_parms_from_quantiles(
  q,
  p = c(0.025, 0.975),
  precision = 0.001,
  derivative.epsilon = 0.001,
  start.with.normal.approx = TRUE,
  start = c(1, 1),
  plot = FALSE
)
```

## Arguments

- q:

  A vector of lower and upper bounds, in which Beta distribution will
  have quantiles specified by `p`. For example, `q=c(0.5,0.99)`

- p:

  The lower and upper quantiles of the range one wants to specify.

- precision:

  Approximation precisions.

- derivative.epsilon:

  Precision of calculating derivative.

- start.with.normal.approx:

  Default is `TRUE`, for normal approximation.

- start:

  Starting values of beta parameters.

- plot:

  Default is `FALSE` to suppress plotting of the beta density,
  otherwise, set to `TRUE`.

## Value

A list containing the selected Beta parameters `a`, and `b`. Other
elements of the list include some details about the computations
involved in finding `a` and `b`.

## Examples

``` r
beta_parms_from_quantiles(c(0.5,0.99))
#> $a
#> [1] 5.967161
#> 
#> $b
#> [1] 1.261544
#> 
#> $last.change
#> [1] -1.493917e-04 -6.358114e-05
#> 
#> $niter
#> [1] 5
#> 
#> $q
#> [1] 0.50 0.99
#> 
#> $p
#> [1] 0.025 0.975
#> 
#> $p.check
#> [1] 0.025 0.975
#> 
```
