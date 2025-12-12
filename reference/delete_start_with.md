# Deletes a pattern from the start of a string, or each of a vector of strings.

`delete_start_with` is used for clean the column names in raw data. For
example, R adds "X" at the start of variable names. This function
deletes "X\_"s from the column names. This can happen if the raw data
have column names such as "`_CASE_ABX`". Check
[`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)
for its actual usage.

## Usage

``` r
delete_start_with(s, vec)
```

## Arguments

- s:

  the pattern (a single string) to be deleted from the start.

- vec:

  a vector of strings with unwanted starting strings (specified by `s`).

## Value

string(s) with deleted patterns from the start.

## Examples

``` r
delete_start_with("X_",c("X_hello"))
#> [1] "hello"
delete_start_with("X_",c("X_hello","hello2"))
#> [1] "hello"  "hello2"
delete_start_with("X_",c("X_hello","hello2","X_hello3"))
#> [1] "hello"  "hello2" "hello3"
```
