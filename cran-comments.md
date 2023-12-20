## Resubmission
This is a resubmission. In this version I have:

* Fixed an issue related to JAGS being updated to 4.3.x 
(Package was archived on CRAN because of a tricky interfacing issue between 
JAGS 4.3.x and R 4.3.x was not corrected in time); I contacted the author of JAGS and he confirmed the issue
was a bug but will not release a debugged version soon. So I went ahead to fix the issue 
in my package; the changes primarily occurred in the `jags2_baker` function which interfaces with
the JAGS. Also JAGS 4.3.x retired a function`cut` so I modified the functoins in `nplcm.R` to accommodate this change.
  - And these fixes passes the suggested checks (see Test `Environments` below). 

* winbuilder gave a note on new submission and possible misspelling - but because
I was intructed to not put single quotes around these names (which would overcome this issue),
I will leave as is.

  - ```
    New submission
    
    Package was archived on CRAN
    
    Possibly misspelled words in DESCRIPTION:
      Deloria (18:9, 19:9)
      Hammitt (18:24)
      Zeger (18:37, 19:28)
  ```

## Reverse dependencies; revdecp_check result `OK: 0` `BROKEN: 0`

This is a precise fix to a previously a short and archived release taken down on 
2022-06-08 (https://cran.r-project.org/web/packages/baker/index.html), 
so there are no reverse dependencies given the short life of the previous release.


## Test environments
* local macOS X Sonoma 14.2 install, R 4.3.2
* ubuntu latest release (on Github-Action), R latest release
* macOS latest release (on Github-Action), R latest release
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note


---

* There are no downstream dependencies.
