# Show function dependencies

Show function dependencies

## Usage

``` r
show_dep(fname, pckg = "package:baker", ...)
```

## Arguments

- fname:

  Character string for one function

- pckg:

  Package name; default is `"package:baker"`

- ...:

  Other parameters accepted by
  [`mvbutils::foodweb()`](https://rdrr.io/pkg/mvbutils/man/foodweb.html)

## Value

A figure showing function dependencies

## Examples

``` r
show_dep("nplcm",ancestor=FALSE)

#> [[1]]
#> plot.foodweb
#> 
#> $x
#> answer
#> 
#> $textcolor
#> textcolor
#> 
#> $boxcolor
#> boxcolor
#> 
#> $xblank
#> xblank
#> 
#> $border
#> border
#> 
#> $color.lines
#> color.lines
#> 
#> $plotmath
#> plotmath
#> 
#> $cex
#> [1] 1
#> 
show_dep("nplcm")

#> [[1]]
#> plot.foodweb
#> 
#> $x
#> answer
#> 
#> $textcolor
#> textcolor
#> 
#> $boxcolor
#> boxcolor
#> 
#> $xblank
#> xblank
#> 
#> $border
#> border
#> 
#> $color.lines
#> color.lines
#> 
#> $plotmath
#> plotmath
#> 
#> $cex
#> [1] 1
#> 
```
