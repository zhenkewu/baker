# Read measurement slices

NB: add example, small data

## Usage

``` r
read_meas_object(object, data)
```

## Arguments

- object:

  Outputs from
  [`make_meas_object()`](https://zhenkewu.com/baker/reference/make_meas_object.md)

- data:

  Raw data with column names `{pathogen name}_{specimen}{test}`

## Value

A list with two elements: `meas`-data frame with measurements;
`position`-see
[`lookup_quality()`](https://zhenkewu.com/baker/reference/lookup_quality.md)

## See also

Other raw data importing functions:
[`extract_data_raw()`](https://zhenkewu.com/baker/reference/extract_data_raw.md)
