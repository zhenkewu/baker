# Create new folder name

Create new folder name

## Usage

``` r
make_foldername(parent_path, parameter_names, parameter_vals, sep = "/")
```

## Arguments

- parent_path:

  The parent directory where to put the new folder

- parameter_names:

  The parameters that distinguish this folder's scenario

- parameter_vals:

  The actual parameter values

- sep:

  file name separator - default to `"/"` for OSX; `"\\"` for Windows.

## Value

A string for folder name

## Examples

``` r
make_foldername("/user",c("theta","alpha","beta"),c(1,2,3))
#> [1] "/user/theta=1_alpha=2_beta=3"

```
