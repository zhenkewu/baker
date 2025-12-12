# get top patterns from a slice of bronze-standard measurement

get top patterns from a slice of bronze-standard measurement

## Usage

``` r
get_top_pattern(BrS_dat, Y, case_status, n_pat, exclude_missing = TRUE)
```

## Arguments

- BrS_dat:

  bronze-standard data, which is usually `data_nplcm$Mobs$MBS[[1]]`

- Y:

  A vector of case/control status: 1 for case; 0 for control

- case_status:

  1 for case; 0 for controls

- n_pat:

  the number of top patterns one wants to show

- exclude_missing:

  DEFAULT is TRUE for excluding any individual with missing
  measurements.

## Value

a list of results: `obs_pat` - observed rates; `pattern_names`;
`exist_other` - if actual no. of patterns is larger than `n_pat`; `N`-
No. of individuals with `Y = case_status`.

## See also

Other exploratory data analysis functions:
[`plot_logORmat()`](https://zhenkewu.com/baker/reference/plot_logORmat.md),
[`show_individual()`](https://zhenkewu.com/baker/reference/show_individual.md),
[`summarize_BrS()`](https://zhenkewu.com/baker/reference/summarize_BrS.md),
[`summarize_SS()`](https://zhenkewu.com/baker/reference/summarize_SS.md),
[`visualize_season()`](https://zhenkewu.com/baker/reference/visualize_season.md)

## Examples

``` r
data(data_nplcm_noreg)
get_top_pattern(data_nplcm_noreg$Mobs$MBS[[1]],data_nplcm_noreg$Y,1,5,FALSE)
#> $obs_pat
#>     100000     001000     101000     110000     010000      other 
#> 0.22666667 0.08666667 0.08666667 0.06666667 0.06333333 0.47000000 
#> 
#> $pattern_names
#> [1] "100000" "001000" "101000" "110000" "010000" "other" 
#> 
#> $exist_other
#> [1] TRUE
#> 
#> $N
#> [1] 300
#> 

data(data_nplcm_noreg)
get_top_pattern(data_nplcm_noreg$Mobs$MBS$MBS1,data_nplcm_noreg$Y,case_status=1,n_pat=5)
#> $obs_pat
#>     100000     001000     101000     110000     010000      other 
#> 0.22666667 0.08666667 0.08666667 0.06666667 0.06333333 0.47000000 
#> 
#> $pattern_names
#> [1] "100000" "001000" "101000" "110000" "010000" "other" 
#> 
#> $exist_other
#> [1] TRUE
#> 
#> $N
#> [1] 300
#> 

```
