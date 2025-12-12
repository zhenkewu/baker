# Convert 0/1 coding to pathogen/combinations

Reverse to [`symb2I()`](https://zhenkewu.com/baker/reference/symb2I.md)

## Usage

``` r
I2symb(binary_code, pathogen_list)
```

## Arguments

- binary_code:

  Binary indicators for pathogens

- pathogen_list:

  The complete list of pathogen names

## Value

The name of pathogen or pathogen combination indicated by "code"

## Examples

``` r
I2symb("001",c("A","B","C"))
#> [1] "C"
I2symb("000",c("A","B","C"))
#> [1] "NoA"
```
