## Resubmission
This is a resubmission. In this version I have:

* updated to fix an issue in the example of `nplcm` (Package was archived on CRAN because of a tricky interfacing issue between JAGS 4.3.x and R 4.3.x)

* [only appeared on Windows Server 2022, R-devel, x64] Lost braces; missing escapes or markup?
  
* There is another note (occasionally) about some doi links being "Service Unavailable". 
  I checked that every link works. They were considered to be slow (by the checking functions) 
  when directed to the page.

* Previously I was asked to ignore the following note about spelling (for author names)
  Possibly misspelled words in DESCRIPTION:
       al (17:49, 19:11)
       et (17:46, 19:8)

## Reverse dependencies

This is a precise fix to a previously a short and archived release taken down on 2022-06-08 (https://cran.r-project.org/web/packages/baker/index.html), so there are no reverse dependencies
given the short life of the previous release.



## Test environments
* local OS X install, R 4.3.1
* ubuntu latest release (on Github-Action), R 4.1.2
* win-builder (devel and release), R 4.1.2

## R CMD check results

0 errors | 0 warnings | 1 note


---

* There are no downstream dependencies.
