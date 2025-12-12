# See if a result folder is obtained by JAGS

See if a result folder is obtained by JAGS

## Usage

``` r
is_jags_folder(DIR_NPLCM)
```

## Arguments

- DIR_NPLCM:

  directory to the folder with results. "mcmc_options.txt" must be in
  the folder.

## Value

TRUE for from JAGS; FALSE otherwise.

## Examples

``` r
is_jags_folder(tempdir()) # just an illustration.
#> [baker] no `mcmc_option.txt` found.
#> [1] FALSE
```
