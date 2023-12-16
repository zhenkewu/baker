## Resubmission
This is a resubmission. In this version I have:

* updated to fix an issue in the example of `nplcm` (Package was archived on CRAN because of a tricky interfacing issue between JAGS 4.3.x and R 4.3.x was not corrected in time)

* [only appeared on Windows Server 2022, R-devel, x64] Lost braces; missing escapes or markup?
  
* There is another note (occasionally) about some doi links being "Service Unavailable". 
  I checked that every link works. They were considered to be slow (by the checking functions) 
  when directed to the page.

* Previously I was asked to ignore the following note about spelling (for author names)
  Possibly misspelled words in DESCRIPTION:
       al (17:49, 19:11)
       et (17:46, 19:8)

## Reverse dependencies; revdecp_check result `OK: 0` `BROKEN: 0`

This is a precise fix to a previously a short and archived release taken down on 2022-06-08 (https://cran.r-project.org/web/packages/baker/index.html), so there are no reverse dependencies
given the short life of the previous release.


## Test environments
* local macOS X Sonoma 14.2 install, R 4.3.2
* ubuntu latest release (on Github-Action), R latest release
* windows latest release (on Github-Action), R latest release
* macOS latest release (on Github-Action), R latest release
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note


---

* There are no downstream dependencies.
