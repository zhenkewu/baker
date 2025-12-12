# Sample a vector of Bernoulli variables.

Sample a vector of Bernoulli variables with higher speed (same length
with `"p"`). The Bernoulli random variables can have different means.

## Usage

``` r
rvbern(p)
```

## Arguments

- p:

  A vector of probabilities, each being the head probability of an
  independent coin toss

## Value

A vector of 1s (head) and 0s (tail)

## Examples

``` r
rvbern(c(0.9,0.1,0.2,0.3))
#> [1] 1 0 1 0
```
