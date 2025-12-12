# Reorder the measurement dimensions to match the order for display

Reorder the measurement dimensions to match the order for display

## Usage

``` r
my_reorder(disp_order, raw_nm)
```

## Arguments

- disp_order:

  The vector of names to be displayed (order matters)

- raw_nm:

  The vector of names from raw measurements (order matters)

## Value

A permuted vector from 1 to `length(raw_nm)`. For example, if its first
element is 3, it means that the 3rd pathogen in `raw_nm` should be
arranged to the first in the raw measurements.

## Examples

``` r
  disp_order <- c("B","E","D","C","F","A")
  raw_nm <- c("C","A","E")
  my_reorder(disp_order,raw_nm)
#> [1] 3 1 2

```
