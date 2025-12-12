# Interpret the specified model structure

`assign_model` translates options specified by a user (e.g., in
`model_options`) into information that can be understood by `baker`.

## Usage

``` r
assign_model(model_options, data_nplcm, silent = TRUE)
```

## Arguments

- model_options:

  See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)
  function.

- data_nplcm:

  Data. See [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md)
  function for data structure.

- silent:

  Default is `TRUE` for no messages; `FALSE` otherwise.

## Value

A list of model specifications:

- `num_slice` A vector counting the No. of measurement slices for each
  level of measurement quality (e.g., MBS, MSS, MGS representing
  Bronze-Standard Measurements - case-control, Silver-Standard
  Measurements and Gold-Standard Measurements - case-only);

- `nested` Local dependence specification for modeling bronze-standard
  data. `TRUE` for nested models (conditional dependence given disease
  class); `FALSE` for non-nested models (conditional independence given
  disease class). One for each BrS slice.

- `regression`

  - `do_reg_Eti` `TRUE` for doing etiology regression. It means let the
    etiology fractions vary with explanatory variables. `FALSE`
    otherwise;

  - `do_reg_FPR` A vector whose names represent the slices of
    bronze-standard data. For each slice of BrS measurements, `TRUE`
    does false positive rate regression. It means the false positive
    rates, estimatable from controls, can vary with covariates; `FALSE`
    otherwise.

  - `is_discrete_predictor` A list of names "Eti", and the names for
    every slice of bronze-standard data. `TRUE` if all predictors are
    discrete; `FALSE` otherwise.

## Details

`assign_model` will be modified to check if data are conformable to
specified model.

## Examples

``` r
cause_list <- c(LETTERS[1:6]) 
J.BrS <- 6
model_options_no_reg <- list(
likelihood   = list(
  cause_list = cause_list,
  k_subclass = 2,
  Eti_formula = ~-1, 
  # no covariate for the etiology regression
  FPR_formula = list(
    MBS1 =   ~-1)    
    # no covariate for the subclass weight regression
),
use_measurements = c("BrS"), 
# use bronze-standard data only for model estimation.
prior= list(
  Eti_prior = overall_uniform(1,cause_list), 
  # Dirichlet(1,...,1) prior for the etiology.
  TPR_prior  = list(BrS = list(
    info  = "informative", # informative prior for TPRs
    input = "match_range", 
    # specify the informative prior for TPRs by specifying a plausible range.
    val = list(MBS1 = list(up =  list(rep(0.99,J.BrS)), 
    # upper ranges: matched to 97.5% quantile of a Beta prior
                           low = list(rep(0.55,J.BrS))))
                           # lower ranges: matched to 2.5% quantile of a Beta prior
  )
  )
)
)     
data("data_nplcm_noreg")

assign_model(model_options_no_reg,data_nplcm_noreg)
#> $num_slice
#> MBS MSS MGS 
#>   1   0   0 
#> 
#> $nested
#> [1] TRUE
#> 
#> $regression
#> $regression$do_reg_Eti
#> [1] FALSE
#> 
#> $regression$do_reg_FPR
#>  MBS1 
#> FALSE 
#> 
#> $regression$is_discrete_predictor
#> $regression$is_discrete_predictor$Eti
#> [1] FALSE
#> 
#> $regression$is_discrete_predictor$FPR
#>  MBS1 
#> FALSE 
#> 
#> 
#> 
#> $BrS_grp
#> [1] FALSE
#> 
#> $SS_grp
#> [1] FALSE
#> 
```
