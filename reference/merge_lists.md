# For a list of many sublists each of which has matrices as its member, we combine across the many sublists to produce a final list

For a list of many sublists each of which has matrices as its member, we
combine across the many sublists to produce a final list

## Usage

``` r
merge_lists(list_of_lists)
```

## Arguments

- list_of_lists:

  a list of sublists

## Value

a list after merge

## See also

Other data operation functions:
[`combine_data_nplcm()`](https://zhenkewu.com/baker/reference/combine_data_nplcm.md),
[`subset_data_nplcm_by_index()`](https://zhenkewu.com/baker/reference/subset_data_nplcm_by_index.md)

## Examples

``` r
DT1 = list(A=1:3,B=letters[1:3])
DT2 = list(A=4:5,B=letters[4:5])
DT3 = list(A=1:4,B=letters[1:4])
DT4 = list(A=4:7,B=letters[4:7])
l = list(DT1,DT2);names(l) <- c("haha","hihi")
l2 = list(DT3,DT4);names(l2) <- c("haha","hihi")
listoflists <- list(l,l2);names(listoflists) <- c("dude1","dude2")
listoflists
#> $dude1
#> $dude1$haha
#> $dude1$haha$A
#> [1] 1 2 3
#> 
#> $dude1$haha$B
#> [1] "a" "b" "c"
#> 
#> 
#> $dude1$hihi
#> $dude1$hihi$A
#> [1] 4 5
#> 
#> $dude1$hihi$B
#> [1] "d" "e"
#> 
#> 
#> 
#> $dude2
#> $dude2$haha
#> $dude2$haha$A
#> [1] 1 2 3 4
#> 
#> $dude2$haha$B
#> [1] "a" "b" "c" "d"
#> 
#> 
#> $dude2$hihi
#> $dude2$hihi$A
#> [1] 4 5 6 7
#> 
#> $dude2$hihi$B
#> [1] "d" "e" "f" "g"
#> 
#> 
#> 
merge_lists(listoflists)
#> $haha
#>      A         B          
#> [1,] integer,3 character,3
#> [2,] integer,4 character,4
#> 
#> $hihi
#>      A         B          
#> [1,] integer,2 character,2
#> [2,] integer,4 character,4
#> 
```
