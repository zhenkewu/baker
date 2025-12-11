## Resubmission

This is a resubmission for the `baker` package (last version 1.0.3) 
that was archived last month because of a check error that was not fixed in time. 
In this version (patch), I have:

* Fixed an issue related to the ggplot2 label requirement change that made
function `plot_check_common_pattern()` produce an error. This removed an Error.
* In the example of the package main function `nplcm()`, removed a 
temporary folder created when executing the example; similarly we fixed the issue
for all other examples that created temporary folders. This removes a Note.

## Reverse dependencies; revdecp_check result `OK: 0` `BROKEN: 0`

## Test environments
* local macOS X Tahoe 26.1 install, R 4.5.2
* ubuntu latest release (on Github-Action), R latest release
* macOS latest release (on Github-Action), R latest release
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

* [a Note only for R Under development (unstable), Windows Server 2022 x64 (build 20348)]
Possibly misspelled words in DESCRIPTION:
  Deloria (18:9, 19:9)
  Hammitt (18:24)
  Zeger (18:37, 19:28)
These words are listed in inst/WORDLIST; other checks do not show this.

---

* There are no downstream dependencies.
