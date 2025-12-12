# combine multiple data_nplcm (useful when simulating data from regression models)

combine multiple data_nplcm (useful when simulating data from regression
models)

## Usage

``` r
combine_data_nplcm(data_nplcm_list)
```

## Arguments

- data_nplcm_list:

  a list of data_nplcm in
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)

## Value

a list with each element resulting from row binding of each
corresponding element in the input `data_nplcm_list`.

## See also

Other data operation functions:
[`merge_lists()`](https://zhenkewu.com/baker/reference/merge_lists.md),
[`subset_data_nplcm_by_index()`](https://zhenkewu.com/baker/reference/subset_data_nplcm_by_index.md)

## Examples

``` r
N=100
Y = rep(c(1,0),times=50) # simulate two cases and two controls.
out_list <- vector("list",length=N)
J = 3                          # number of causes
cause_list = c(LETTERS[1:J])   # cause list
K = 2                          # number of subclasses
lambda = c(.8,.2)                # subclass weights for control group
eta = c(.9,.1)                   # subclass weights for case group

for (i in 1:N){
  #setup parameters for the present individual:
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
    Nd      =     1,
    Nu      =     1
  )
  simu_out   <- simulate_nplcm(set_parameter)
  out <- simu_out$data_nplcm
  out_list[[i]] <- out
}

# extract cases and controls and combine all the data into one:
data_nplcm_list <- lapply(1:N, function(s) subset_data_nplcm_by_index(out_list[[s]],2-Y[s]))
data_nplcm_unordered      <- combine_data_nplcm(data_nplcm_list)
```
