## Resubmission
This is a resubmission. In this version I have:

* fixed Author@R issue; removed all functions based on a code snippet online (related to author issue; now all contributions/copyright holders are from the currently listed authors.); in particular, the function at issue`sort_data_frame()` is now completely removed.

* removed "An R package for" in tht title

* removed the quotes around et and al

* added \value in all the 14 .Rd files that previously did not have this tag;
(manually checked using `grep -L -F "\value{" man/*.Rd`, which did not return any files).

* fixed issues related to examples in unexported functions, I kept the functions unexported and removed examples.

* replaced \dontrun with \donttest for examples that involve
`>`5 seconds because of posterior MCMC sampling steps; for other previous instances of \dontrun
we have unwrapped them by providing fast examples

* checked functions; did not find any actions that would write in user's home filespace

* fixed the issue related to `par()`; now all changes by plotting functions in examples
and vignettes because of`par()` are reset after the plots are rendered

* removed all instances `rm(list=ls())` in examples and vignettes

* removed all instances of `installed.packages()`


## Test environments
* local OS X install, R 4.0.2
* ubuntu latest release (on Github-Action), R 4.1.2
* win-builder (devel and release), R 4.1.2

## R CMD check results

0 errors | 0 warnings | 1 note

* New submission

* [only appeared on Windows Server 2022, R-devel, 64 bit]checking for detritus in the temp directory ... NOTE
Found the following files/directories:
  'lastMiKTeXException'
  We checked online; it seems to be platform-specific and related to files produced during testing; we did not find any place in the source files where interactive use of the vignette was invoked. This note did not appear in other platforms.

## Reverse dependencies

This is a new release, so there are no reverse dependencies.

---

* There are no downstream dependencies.
