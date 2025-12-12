# get an individual's data from the output of [`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

get an individual's data from the output of
[`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

## Usage

``` r
show_individual(data_nplcm, ID)
```

## Arguments

- data_nplcm:

  data for fitting nplcm; See
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- ID:

  patient id: `patid`.

## Value

a list with the inquired patient's data

## See also

Other exploratory data analysis functions:
[`get_top_pattern()`](https://zhenkewu.com/baker/reference/get_top_pattern.md),
[`plot_logORmat()`](https://zhenkewu.com/baker/reference/plot_logORmat.md),
[`summarize_BrS()`](https://zhenkewu.com/baker/reference/summarize_BrS.md),
[`summarize_SS()`](https://zhenkewu.com/baker/reference/summarize_SS.md),
[`visualize_season()`](https://zhenkewu.com/baker/reference/visualize_season.md)

## Examples

``` r
data(data_nplcm_noreg)
data_nplcm_noreg$X$patid <- paste("PAT",1:length(data_nplcm_noreg$Y0),sep="")
data_nplcm_noreg$X <- as.data.frame(data_nplcm_noreg$X)
subset_data_nplcm_by_index(data_nplcm_noreg,which(data_nplcm_noreg$X$patid%in%c("PAT12","PAT408")))
#> $Mobs
#> $Mobs$MBS
#> $Mobs$MBS$MBS1
#> [1] A B C D E F
#> <0 rows> (or 0-length row.names)
#> 
#> 
#> $Mobs$MSS
#> $Mobs$MSS$MSS1
#> [1] A B
#> <0 rows> (or 0-length row.names)
#> 
#> 
#> 
#> $Y
#> numeric(0)
#> 
#> $X
#> [1] patid
#> <0 rows> (or 0-length row.names)
#> 
data_nplcm_noreg$X <- NULL
```
