# Takes any number of R objects as arguments and returns a list whose names are derived from the names of the R objects.

Roger Peng's listlabeling challenge from
<https://simplystatistics.tumblr.com/post/11988685443/computing-on-the-language>.
Code copied from
<https://gist.github.com/ajdamico/1329117/0134148987859856fcecbe4446cfd37e500e4272>

## Usage

``` r
make_list(...)
```

## Arguments

- ...:

  any R objects

## Value

a list as described above

## Examples

``` r
#create three example variables for a list
x <- 1
y <- 2
z <- "hello"
#display the results
make_list( x , y , z )
#> $x
#> [1] 1
#> 
#> $y
#> [1] 2
#> 
#> $z
#> [1] "hello"
#> 
```
