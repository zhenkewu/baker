# compute positive rates for nested model with subclass mixing weights that are the same across `Jcause` classes for each person (people may have different weights.)

The array version of this function
([compute_marg_PR_nested_reg_array](https://zhenkewu.com/baker/reference/compute_marg_PR_nested_reg_array.md))
is used in
[plot_etiology_regression](https://zhenkewu.com/baker/reference/plot_etiology_regression.md)

## Usage

``` r
compute_marg_PR_nested_reg(ThetaBS, PsiBS, pEti_mat, subwt_mat, case, template)
```

## Arguments

- ThetaBS:

  True positive rates for `JBrS` measures (rows) among `K` subclasses
  (columns)

- PsiBS:

  False positive rates; dimension same as above

- pEti_mat:

  a matrix of etiology pies for `N` subjects (rows) and `Jcause` causes
  (columns) rows sum to ones.

- subwt_mat:

  a matrix of subclass weights for cases and controls. `N` by `K`. Rows
  sum to ones.

- case:

  a N-vector of `1`s (cases) and `0`s (controls)

- template:

  a binary matrix with `Jcause+1` rows (`Jcause` classes of cases and
  `1` class of controls) and `JBrS` columns for the Bronze-standard
  measurement (say, pick one type/slice). The ones in each row indicate
  the measurements that will show up more frequently in cases given the
  cause.

## Value

a matrix of values between `0` and `1` (need not to have row sums of
ones); of dimension (number of subjects, dimension of the
bronze-standard measurement slice).
