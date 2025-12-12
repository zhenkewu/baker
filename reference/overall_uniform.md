# specify overall uniform (symmetric Dirichlet distribution) for etiology prior

specify overall uniform (symmetric Dirichlet distribution) for etiology
prior

## Usage

``` r
overall_uniform(alpha, cause_list)
```

## Arguments

- alpha:

  any positive number, usually 1.

- cause_list:

  a list of latent status

## Value

a vector of length `length(cause_list)`

## See also

Other prior specification functions:
[`set_prior_tpr_BrS_NoNest()`](https://zhenkewu.com/baker/reference/set_prior_tpr_BrS_NoNest.md),
[`set_prior_tpr_SS()`](https://zhenkewu.com/baker/reference/set_prior_tpr_SS.md)

## Examples

``` r
overall_uniform(1,c("A","B","C"))
#> [1] 1 1 1
```
