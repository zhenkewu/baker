# compute positive rates for nested model with subclass mixing weights that are the same across `Jcause` classes for each person (people may have different weights.)

This is an array-version of
[compute_marg_PR_nested_reg](https://zhenkewu.com/baker/reference/compute_marg_PR_nested_reg.md).
This is used in
[plot_etiology_regression](https://zhenkewu.com/baker/reference/plot_etiology_regression.md)

## Usage

``` r
compute_marg_PR_nested_reg_array(
  ThetaBS_array,
  PsiBS_array,
  pEti_mat_array,
  subwt_mat_array,
  case,
  template
)
```

## Arguments

- ThetaBS_array:

  An array of: True positive rates for JBrS measures (rows) among K
  subclasses (columns)

- PsiBS_array:

  An array of: False positive rates; dimension same as above

- pEti_mat_array:

  An array of: a matrix of etiology pies for N subjects (rows) and
  Jcause causes (columns) rows sum to ones.

- subwt_mat_array:

  An array of: a matrix of subclass weights for cases and controls. N
  by K. Rows sum to ones.

- case:

  a N-vector of 1s (cases) and 0s (controls)

- template:

  a binary matrix with Jcause+1 rows (Jcause classes of cases and 1
  class of controls) and JBrS columns for the Bronze-standard
  measurement (say, pick one type/slice). The ones in each row indicate
  the measurements that will show up more frequently in cases given the
  cause.

## Value

An array of: a matrix of values between 0 and 1 (need not to have row
sums of ones); of dimension (number of subjects, dimension of the
bronze-standard measurement slice).
