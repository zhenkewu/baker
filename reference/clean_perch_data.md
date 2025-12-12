# Clean PERCH data

`clean_perch_data` transforms a raw data table (row for subjects, column
for variables - named as `{pathogen name}_{specimen}{test}` for lab
tests or other covariates) into a list. It is designed for PERCH data
format.

## Usage

``` r
clean_perch_data(clean_options)
```

## Arguments

- clean_options:

  The list of options for cleaning PERCH data. Its elements are defined
  as follows:

  `raw_meas_dir`

  :   : The file path to the raw data;

  `case_def`

  :   : Variable name in raw data for **case** definition;

  `case_def_val`

  :   : The value for **case** definition;

  `ctrl_def`

  :   : Variable name in raw data for **control** definition;

  `ctrl_def_val`

  :   : The value for **control** definition;

  `X_strat`

  :   : A vector of variable names for stratifying the data to perform
      SEPARATE analyses;

  `X_strat_val`

  :   : A list of values for `X_strat`. The output data only have
      individuals with `identical(X_strat,X_strat_val)==TRUE`. To
      perform analysis on a single site, say `"02GAM"`, use
      `X_strat="newSITE"` and `X_strat_val=list("02GAM")`;

  `BrS_objects`

  :   : A list of BrS objects built by
      [`make_meas_object()`](https://zhenkewu.com/baker/reference/make_meas_object.md);

  `SS_objects`

  :   : A list of SS objects built by
      [`make_meas_object()`](https://zhenkewu.com/baker/reference/make_meas_object.md);

  `X_extra`

  :   : A vector of covariate names for regression or visualization;

  `patho_taxo_dir`

  :   : The file path to the pathogen category or taxonomy information
      (.csv). The information should be as complete as possible for a
      particular analysis. If not, the pathogen without taxonomy
      information could not be assigned to bacterial or viral groups
      (see `plot_group_etiology()`); See `assign_taxo_cause_list()` that
      requires this taxonomy information.

## Value

A List: `list(Mobs,Y,X)`

- `Mobs` A list of bronze- (`MBS`), silver- (`MSS`), and gold-standard
  (`MGS`, if available) measurements. See the formats of these
  measurements in
  [`extract_data_raw()`](https://zhenkewu.com/baker/reference/extract_data_raw.md).

- `Y` 1 for case; 0 for control;

- `X` Data frame of covariates for cases and controls. The covariate
  names are specified in `X_extra`;

This function does not re-order pathogens that only have silver-standard
data.

## See also

[make_meas_object](https://zhenkewu.com/baker/reference/make_meas_object.md)
for wrapping information about a particular type of measurement;
[extract_data_raw](https://zhenkewu.com/baker/reference/extract_data_raw.md)
for reading raw data table and organizing them into `data_nplcm` format.
Also see
[clean_combine_subsites](https://zhenkewu.com/baker/reference/clean_combine_subsites.md)
for combining subsites and
[parse_date_time](https://lubridate.tidyverse.org/reference/parse_date_time.html)
for parsing date.
