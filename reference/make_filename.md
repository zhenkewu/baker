# Create new file name

Create new file name

## Usage

``` r
make_filename(parameter_names, parameter_vals, format)
```

## Arguments

- parameter_names:

  The parameters that distinguish this folder's scenario

- parameter_vals:

  The actual parameter values

- format:

  The suffix ".XXX" in the end to specify the file format

## Value

A string for file name

## Examples

``` r
make_filename(c("theta","alpha"),c(0.9,2),"csv")
#> [1] "theta=0.9_alpha=2.csv"
```
