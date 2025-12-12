# baker (development version)

# baker 1.0.4

* Fixed an issue related to the ggplot2 label requirement change that made
function `plot_check_common_pattern()` produce an error. This removed an Error.
* In the example of the package main function `nplcm()`, removed a 
temporary folder created when executing the example; similarly we fixed the issue
for all other examples that created temporary folders. This removes a Note.

# baker 1.0.3

Fixed an issue related to `is.R` deprecation for R 4.4.0 and above.

# baker 1.0.2

* fix an issue related to JAGS 4.3.x removed `cut()` function and only recognizes
`.Dim` not `dim` in array data objects stored in `jagsdata.txt` now. 

# baker 1.0.0

* first public CRAN release



