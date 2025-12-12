# test if a formula has terms not created by \[s_date_Eti() or [`s_date_FPR()`](https://zhenkewu.com/baker/reference/s_date_FPR.md)

test if a formula has terms not created by \[s_date_Eti() or
[`s_date_FPR()`](https://zhenkewu.com/baker/reference/s_date_FPR.md)

## Usage

``` r
has_non_basis(form)
```

## Arguments

- form:

  a formula

## Value

logical `TRUE` (if having terms not created by \[s_date_Eti() or
[`s_date_FPR()`](https://zhenkewu.com/baker/reference/s_date_FPR.md));
`FALSE` otherwise.

## Examples

``` r
form1 <- as.formula(~ -1+s_date_FPR(DATE,Y,basis = "ps",10) + as.factor(SITE))
form2 <- as.formula(~ -1+s_date_FPR(DATE,Y,basis = "ps",10))
form3 <- as.formula(~ s_date_FPR(DATE,Y,basis = "ps",10))

has_non_basis(form1)
#> [1] TRUE
has_non_basis(form2)
#> [1] FALSE
has_non_basis(form3)
#> [1] TRUE
```
