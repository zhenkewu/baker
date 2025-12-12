# Convert names of pathogen/combinations into 0/1 coding

Convert names of pathogen/combinations into 0/1 coding

## Usage

``` r
symb2I(pathogen_name, pathogen_list)
```

## Arguments

- pathogen_name:

  The allowed pathogen name (can be a combination of pathogens in
  "pathlist")

- pathogen_list:

  The complete list of pathogen names

## Value

A 1 by length(pathlist) matrix of binary code (usually for pathogen
presence/absence)

## Examples

``` r
symb2I("A",c("A","B","C"))
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
symb2I("A+B",c("A","B","C"))
#>      [,1] [,2] [,3]
#> [1,]    1    1    0
symb2I("NoA",c("A","B","C"))
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
symb2I(c("A","B+C"),c("A","B","C")) # gives a 2 by 3 matrix.
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    1

```
