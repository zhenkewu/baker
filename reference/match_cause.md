# Match latent causes that might have the same combo but different specifications

@details In our cause_list, "A+B" represents the same cause as "B+A". It
is used for plotting side-by-side posterior sample comparisons

## Usage

``` r
match_cause(pattern, vec)
```

## Arguments

- pattern:

  a vector of latent cause names, e.g., from a particular fit

- vec:

  a vector of latent cause names, e.g., usually a union of cause names
  from several model fits. Usually, it is also the display order that
  one wants to show.

## Value

A vector of length `length(vec)`; `NA` means no pattern matches vec; 1
at position 10 means the first element of `pattern` matches the 10th
element of `vec`.

## Examples

``` r
pattern <- c("X+Y","A+Z","C")
vec     <- c(LETTERS[1:26],"Y+Z","Y+X","Z+A")
match_cause(pattern,vec)
#>  [1] NA NA  3 NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> [26] NA NA  1  2
```
