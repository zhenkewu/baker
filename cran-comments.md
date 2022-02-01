## Resubmission
This is a resubmission. In this version I have:

* exported functions that used `baker:::`

* removed all occurrences of `rm(list = ls())`


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
