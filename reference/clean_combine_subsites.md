# Combine subsites in raw PERCH data set

In the Actual PERCH data set, a study site may have multiple subsites.
`clean_combine_subsites` combines all the study subjects from the same
site.

## Usage

``` r
clean_combine_subsites(raw_meas_dir, subsites_list, newsites_vec)
```

## Arguments

- raw_meas_dir:

  The file path to the raw data file (.csv)

- subsites_list:

  The list of subsite group names. Each group is a vector of subsites to
  be combined

- newsites_vec:

  A vector of new site names. It has the same length as
  `"subsites_list"`

## Value

A data frame with combined sites
