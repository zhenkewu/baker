# check if a list has elements all of length one

check if a list has elements all of length one

## Usage

``` r
is_length_all_one(x)
```

## Arguments

- x:

  a list

## Value

TRUE or FALSE

## Examples

``` r
l = list(a = 5, b = 1:2)
is_length_all_one(l) # FALSE
#> [1] FALSE
l = list(a = 5, b = 1)
is_length_all_one(l) # TRUE
#> [1] TRUE
```
