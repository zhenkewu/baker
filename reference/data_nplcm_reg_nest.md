# Simulated dataset that is structured in the format necessary for an [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md) with regression

Data set for illustrating regression functionalities

## Usage

``` r
data("data_nplcm_reg_nest")
```

## Format

A list containing three items

- Mobs:

  BrS level measurements: N = 1,200 (half cases and half controls); one
  slice of BrS measurements (6 dimensional, A-F); one slice of SS
  measurements (2 dimensional, A and B)

- Y:

  case-control status

- X:

  matrix of covariates (N by 4); columns: SITE (1 and 2, each with 600
  subjects), DATE (index from 1:300), std_date (standardized DATE),
  ENRLDATE (actual dates)

## Value

No returned value; just loading data into the working space.
