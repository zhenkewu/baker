# check existence and create folder if non-existent

check existence and create folder if non-existent

## Usage

``` r
check_dir_create(path)
```

## Arguments

- path:

  Folder path to check and create if not there.

## Value

the same returned values for
[`dir.create()`](https://rdrr.io/r/base/files2.html)

## Examples

``` r
# \donttest{
run_example <- function(){
 xx <- tempdir()
 check_dir_create(xx)
 # on.exit(unlink(xx, recursive = TRUE), add = TRUE) 
}
run_example()
#> NULL
# }
```
