# calculate pairwise log odds ratios

Case at upper triangle; control at lower triangle

## Usage

``` r
logOR(MBS.case, MBS.ctrl)
```

## Arguments

- MBS.case:

  Case Bronze-Standard (BrS) data; rows for case subjects; columns
  contain `JBrS` measurements

- MBS.ctrl:

  Control Bronze-Standard (BrS) data; rows for control subjects; columns
  contain `JBrS` measurements

## Value

a list of two elements: `logOR` (`JBrS` by `JBrS` matrix of log odds
ratios for each pair among `JBrS` measurements) and `logOR.se` ( same
dimension as `logOR`, but representing the standard errors of the
corresponding estimated log odds ratios in `logOR`).
