# Make measurement slice

Wrap the information about a particular type of measurement, e.g.,
NPPCR. NB: add example! copy some from the vignette file.

## Usage

``` r
make_meas_object(patho, specimen, test, quality, cause_list, sep_char = "_")
```

## Arguments

- patho:

  A vector of pathogen names

- specimen:

  Specimen name

- test:

  Test name

- quality:

  Quality category: any of "BrS", "SS" or "GS".

- cause_list:

  The vector of potential latent status

- sep_char:

  a character string that separate the pathogen names and the
  specimen-test pair; Default to `"_"`

## Value

A list with measurement information

- `quality` same as argument

- `patho` same as argument

- `name_in_data` the names used in the raw data to locate these
  measurements

- `template` a mapping from `patho` to `cause_list`.
  `NROW = length(cause_list)+1`; `NCOL = length(patho)`. This value is
  crucial in model fitting to determine which measurements are
  informative of a particular category of latent status.

- `specimen` same as argument

- `test` same as argument

- `nm_spec_test` paste `specimen` and `test` together

## See also

[`make_template()`](https://zhenkewu.com/baker/reference/make_template.md)

## Examples

``` r
make_meas_object(
patho = c("A","B","C","D","E","F"), 
specimen = "MBS",
test = "1",
quality = "BrS", 
cause_list = c("A","B","C","D","E"))
#> $quality
#> [1] "BrS"
#> 
#> $patho
#> [1] "A" "B" "C" "D" "E" "F"
#> 
#> $name_in_data
#> [1] "A_MBS1" "B_MBS1" "C_MBS1" "D_MBS1" "E_MBS1" "F_MBS1"
#> 
#> $template
#>      [,1] [,2] [,3] [,4] [,5] [,6]
#> [1,]    1    0    0    0    0    0
#> [2,]    0    1    0    0    0    0
#> [3,]    0    0    1    0    0    0
#> [4,]    0    0    0    1    0    0
#> [5,]    0    0    0    0    1    0
#> [6,]    0    0    0    0    0    0
#> 
#> $specimen
#> [1] "MBS"
#> 
#> $test
#> [1] "1"
#> 
#> $nm_spec_test
#> [1] "MBS1"
#> 
```
