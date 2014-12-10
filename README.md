* How to install?
```r
install.packages(devtools)
devtools::install_github("zhenkewu/nplcm")
```

* Why should someone use your package?

* How does it compare to other existing solutions?

* What are the main functions?

* Platform:
  * currently implemented for Windows system, 7, or 8, because `nplcm` uses WinBUGS
  to fit the models;
  * for Mac OS X system, one can install the package and study the components of
  each function. If you find package `ks` cannot be loaded due to failure of 
  loading package `rgl`, follow the following steps:
    * install X11 by going [here](http://xquartz.macosforge.org/trac/wiki/X112.7.7)
    * `install.packages("http://cran.r-project.org/src/contrib/rgl_0.95.1158.tar.gz"    ,repo=NULL,type="source")`

