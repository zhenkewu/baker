# get individual data

get individual data

## Usage

``` r
get_individual_data(i, data_nplcm)
```

## Arguments

- i:

  index of individual as appeared in `data_nplcm`

- data_nplcm:

  the data for nplcm; see
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

## Value

a list of the same structure as `data_nplcm`; just with one row of
values

## Examples

``` r
data(data_nplcm_noreg)
get_individual_data(2,data_nplcm_noreg)
#> $Mobs
#> $Mobs$MBS
#> $Mobs$MBS$MBS1
#>   A B C D E F
#> 2 1 0 1 0 0 0
#> 
#> 
#> $Mobs$MSS
#> $Mobs$MSS$MSS1
#>   A B
#> 2 0 0
#> 
#> 
#> $Mobs$MGS
#> NULL
#> 
#> 
#> $X
#> NULL
#> 
#> $Y
#> [1] 1
#> 
```
