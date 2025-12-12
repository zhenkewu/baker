# summarize bronze-standard data

summarize bronze-standard data

## Usage

``` r
summarize_BrS(BrS_dat, Y)
```

## Arguments

- BrS_dat:

  bronze-standard data, which is usually `data_nplcm$Mobs$MBS[[1]]`

- Y:

  A vector of case/control status: 1 for case; 0 for control

## Value

a list of summaries for BrS data

## See also

Other exploratory data analysis functions:
[`get_top_pattern()`](https://zhenkewu.com/baker/reference/get_top_pattern.md),
[`plot_logORmat()`](https://zhenkewu.com/baker/reference/plot_logORmat.md),
[`show_individual()`](https://zhenkewu.com/baker/reference/show_individual.md),
[`summarize_SS()`](https://zhenkewu.com/baker/reference/summarize_SS.md),
[`visualize_season()`](https://zhenkewu.com/baker/reference/visualize_season.md)

## Examples

``` r
data(data_nplcm_noreg)
summarize_BrS(data_nplcm_noreg$Mobs$MBS[[1]], data_nplcm_noreg$Y)
#> $Nd
#> [1] 300
#> 
#> $Nu
#> [1] 300
#> 
#> $MBS_mean
#>         [,1]    [,2]    [,3]    [,4]    [,5]    [,6]
#> case    0.60 0.31333 0.37000 0.16333 0.11667 0.13667
#> control 0.21 0.21667 0.23667 0.12000 0.13333 0.13000
#> 
#> $MBS_q1
#>              [,1]      [,2]      [,3]       [,4]       [,5]       [,6]
#> case    0.5436223 0.2634318 0.3173072 0.12559292 0.08477613 0.10209823
#> control 0.1675776 0.1736517 0.1919746 0.08764294 0.09919185 0.09629289
#> 
#> $MBS_q2
#>              [,1]     [,2]      [,3]      [,4]      [,5]      [,6]
#> case    0.6538491 0.367955 0.4259800 0.2095867 0.1582502 0.1804223
#> control 0.2597553 0.266846 0.2880173 0.1619657 0.1767463 0.1730629
#> 
#> $MBS_nm
#> [1] "A" "B" "C" "D" "E" "F"
#> 
```
