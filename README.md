A R Package for Fitting Bayesian [Nested Partially Latent Class Models](http://biostats.bepress.com/jhubiostat/paper276/) 

[![Build Status](https://travis-ci.org/zhenkewu/nplcm.svg?branch=master)](https://travis-ci.org/zhenkewu/nplcm)

How to install?
--------------
```r
install.packages("devtools",repos="http://watson.nci.nih.gov/cran_mirror/")
devtools::install_github("zhenkewu/nplcm")
```

Why should someone use `nplcm`?
-------------------------------------

- To study disease etiology from case-control data from multiple sources that have measurement errors. If you are interested in estimating the population etiology pie (fraction), and the probability of each cause for individual case, try `nplcm`.

- Reference publication can be found [here](http://onlinelibrary.wiley.com/doi/10.1111/rssc.12101/abstract) and [here](http://biostats.bepress.com/jhubiostat/paper276/).

How does it compare to other existing solutions?
------------------------------------------------
- Acknowledges various levels of measurement errors and combines multiple sources
of data.

What are the main functions?
-----------------------------
- `nplcm` that fits the model with or without covariates.

Platform
---------
- Currently implemented for Windows system, 7 or 8, because `nplcm` uses WinBUGS
  to fit the models;
  
- The `.bug` model files are included [here](https://github.com/zhenkewu/bugs.models/tree/master/nplcm);

- For Mac OS X system, one can install the package and study the components of
  each function. If you find package `ks` cannot be loaded due to failure of 
  loading package `rgl`, follow the following steps:
  
    1. install X11 by going [here](http://xquartz.macosforge.org/trac/wiki/X112.7.7);
    
    2. `install.packages("http://download.r-forge.r-project.org/src/contrib/rgl_0.95.1200.tar.gz",repo=NULL,type="source")`

