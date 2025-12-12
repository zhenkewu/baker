# get unique causes, regardless of the actual order in combo

get unique causes, regardless of the actual order in combo

## Usage

``` r
unique_cause(cause_vec)
```

## Arguments

- cause_vec:

  a vector of characters with potential combo repetitions written in
  scrambled orders separated by "+"

## Value

a vector of characters with unique meanings for latent causes

## Examples

``` r
x <- c("A","B","A","CC+DD","DD+CC","E+F+G","B")
unique_cause(x)
#> [1] "A"     "B"     "CC+DD" "E+F+G"
```
