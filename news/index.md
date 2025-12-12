# Changelog

## baker (development version)

## baker 1.0.4

CRAN release: 2025-12-11

- Fixed an issue related to the ggplot2 label requirement change that
  made function
  [`plot_check_common_pattern()`](https://zhenkewu.com/baker/reference/plot_check_common_pattern.md)
  produce an error. This removed an Error.
- In the example of the package main function
  [`nplcm()`](https://zhenkewu.com/baker/reference/nplcm.md), removed a
  temporary folder created when executing the example; similarly we
  fixed the issue for all other examples that created temporary folders.
  This removes a Note.

## baker 1.0.3

CRAN release: 2024-01-30

Fixed an issue related to `is.R` deprecation for R 4.4.0 and above.

## baker 1.0.2

CRAN release: 2023-12-21

- fix an issue related to JAGS 4.3.x removed
  [`cut()`](https://rdrr.io/r/base/cut.html) function and only
  recognizes `.Dim` not `dim` in array data objects stored in
  `jagsdata.txt` now.

## baker 1.0.0

CRAN release: 2022-02-02

- first public CRAN release
