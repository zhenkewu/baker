# subset data from the output of [`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

It is particularly useful in simulating data from a regression model
where one generates a case and control at a particular covariate value,
and just choose a case or control to retain in the simulated data.

## Usage

``` r
subset_data_nplcm_by_index(data_nplcm, index)
```

## Arguments

- data_nplcm:

  data for fitting nplcm; See
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

- index:

  a vector of indices indicating the observations you hope to subset; it
  will subset in all the sublists of data_nplcm

## Value

a list with the requested data, in the order determined by 'index'

## See also

Other data operation functions:
[`combine_data_nplcm()`](https://zhenkewu.com/baker/reference/combine_data_nplcm.md),
[`merge_lists()`](https://zhenkewu.com/baker/reference/merge_lists.md)

## Examples

``` r
J = 3                          # number of causes
cause_list = c(LETTERS[1:J])   # cause list
K = 2                          # number of subclasses
lambda = c(1,0)                # subclass weights for control group
eta = c(1,0)                   # subclass weights for case group

# setup parameters for the present individual:
set_parameter <- list(
 cause_list      = cause_list,
 etiology        = c(0.5,0.2,0.3), # only meaningful for cases 
 pathogen_BrS    = LETTERS[1:J],
 pathogen_SS     = LETTERS[1:2],
 meas_nm         = list(MBS = c("MBS1"),MSS=c("MSS1")),
 Lambda          = lambda,         # for BrS   
 Eta             = t(replicate(J,eta)),  # case subclass weight for BrS
 PsiBS           = cbind(c(0.15,0.3,0.35),   
                         c(0.25,0.2,0.15)), # FPR
 PsiSS           = cbind(rep(0,J),rep(0,J)),
 ThetaBS         = cbind(c(0.95,0.9,0.85),    # TPR
                         c(0.95,0.9,0.85)),
 ThetaSS         = cbind(c(0.25,0.10),
                         c(0.25,0.10)),
 Nd      =     5,
 Nu      =     3 
)
simu_out   <- simulate_nplcm(set_parameter)
out <- simu_out$data_nplcm
out
#> $Mobs
#> $Mobs$MBS
#> $Mobs$MBS$MBS1
#>   A B C
#> 1 1 1 0
#> 2 0 1 0
#> 3 1 0 1
#> 4 1 0 0
#> 5 0 1 1
#> 6 0 0 1
#> 7 0 0 0
#> 8 0 0 0
#> 
#> 
#> $Mobs$MSS
#> $Mobs$MSS$MSS1
#>   A B
#> 1 0 0
#> 2 0 0
#> 3 0 0
#> 4 0 0
#> 5 0 0
#> 6 0 0
#> 7 0 0
#> 8 0 0
#> 
#> 
#> $Mobs$MGS
#> NULL
#> 
#> 
#> $Y
#> [1] 1 1 1 1 1 0 0 0
#> 
#> $X
#> NULL
#> 
subset_data_nplcm_by_index(out,c(1,4,5))
#> $Mobs
#> $Mobs$MBS
#> $Mobs$MBS$MBS1
#>   A B C
#> 1 1 1 0
#> 4 1 0 0
#> 5 0 1 1
#> 
#> 
#> $Mobs$MSS
#> $Mobs$MSS$MSS1
#>   A B
#> 1 0 0
#> 4 0 0
#> 5 0 0
#> 
#> 
#> 
#> $Y
#> [1] 1 1 1
#> 
#> $X
#> NULL
#> 
subset_data_nplcm_by_index(out,2)
#> $Mobs
#> $Mobs$MBS
#> $Mobs$MBS$MBS1
#>   A B C
#> 2 0 1 0
#> 
#> 
#> $Mobs$MSS
#> $Mobs$MSS$MSS1
#>   A B
#> 2 0 0
#> 
#> 
#> 
#> $Y
#> [1] 1
#> 
#> $X
#> NULL
#> 
```
