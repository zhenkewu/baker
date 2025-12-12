# Import Raw PERCH Data `extract_data_raw` imports and converts the raw data to analyzable format

Import Raw PERCH Data

`extract_data_raw` imports and converts the raw data to analyzable
format

## Usage

``` r
extract_data_raw(
  dat_prepared,
  strat_nm,
  strat_val,
  meas_object,
  extra_covariates = NULL
)
```

## Arguments

- dat_prepared:

  The data set prepared in `clean_perch_data`.

- strat_nm:

  The vector of covariate names to separately extract data. For example,
  in PERCH data cleaning, `X = c("newSITE","CASECONT")`.

- strat_val:

  The list of covariate values to stratify data. Each element
  corresponds to elements in `X`. For example, in PERCH data cleaning,
  `Xval = list("02GAM","1")`.

- meas_object:

  A list of bronze-standard or silver-standard measurement objects made
  by function
  [`make_meas_object()`](https://zhenkewu.com/baker/reference/make_meas_object.md).

- extra_covariates:

  The vector of covariate name for regression purposes. The default is
  NULL, which means not reading in any covariate.

## Value

A list of data.

- Mobs:

  MBS

  :   A list of Bronze-Standard (BrS) measurements. The names of the
      list take the form of `specimen`\_`test`. Each element of the list
      is a data frame. The rows of the data frame are for subjects; the
      columns are for measured pathogens.

  MSS

  :   A list of Silver-Standard (SS) measurements. The formats are the
      same as `MBS` above.

  MGS

  :   A list of Gold-Standard (GS) measurements. It equals `NULL` if no
      GS data exist.

- X:

  A data frame with columns specified by `extra_covariates`.

## See also

[`clean_perch_data()`](https://zhenkewu.com/baker/reference/clean_perch_data.md)

Other raw data importing functions:
[`read_meas_object()`](https://zhenkewu.com/baker/reference/read_meas_object.md)
