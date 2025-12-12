# Check if covariates are discrete

`is_discrete` checks if the specified covariates could be regarded as
discrete variables.

## Usage

``` r
is_discrete(X, X_reg)
```

## Arguments

- X:

  A data frame of covariates

- X_reg:

  The vector of covariates that will stratify the analyses. These
  variables have to be categorical. Or a formula (can be tested by
  `is.formula` in `plyr`), e.g.,
  `~as.factor(SITE8) + as.factor(AGECAT > 1)`.

## Value

`TRUE` for all being discrete; `FALSE` otherwise.

## Details

Note that this function should be used with caution. It used
\$\$nrow(X)/nrow(unique(X\[,X_reg,drop=FALSE\]))\>10\$\$ as an *ad hoc*
criterion. It is not the same as `is.discrete()` in `plyr`
