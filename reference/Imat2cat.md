# Convert a matrix of binary indicators to categorical variables

Convert a matrix of binary indicators to categorical variables

## Usage

``` r
Imat2cat(binary_mat, cause_list, pathogen_list)
```

## Arguments

- binary_mat:

  The matrix of binary indicators. Rows for subjects, columns for
  pathogens in the `"pathogen.list"`

- cause_list:

  The list of causes

- pathogen_list:

  The complete list of pathogen names

## Value

A vector of categorical variables. Its length equals the length of
`"allowed.list"`

## Examples

``` r
Imat2cat(rbind(diag(3),c(1,1,0),c(0,0,0)),c("A","B","C","A+B","NoA"),c("A","B","C"))
#> 100 010 001 110 000 
#>   1   2   3   4   5 
```
