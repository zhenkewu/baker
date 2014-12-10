How to install?
--------------
```r
install.packages(devtools)
devtools::install_github("zhenkewu/nplcm")
```

Why should someone use your package?
-------------------------------------

- To study disease etiology from case-control data from multiple sources and have measurement errors. If you are interested in estimating population etiology pie (fraction), and the probability of each cause for individual case, try `nplcm`.

- Reference can be found [here](http://arxiv.org/abs/1411.5774).

How does it compare to other existing solutions?
------------------------------------------------
- Accknowledges various levels of measurement errors and combines multiple sources
of data.

What are the main functions?
-----------------------------
- `nplcm` that fits the model with or without covarites.

Platform
---------
- Currently implemented for Windows system, 7 or 8, because `nplcm` uses WinBUGS
  to fit the models;
  
- The `.bug` model files are included [here](https://github.com/zhenkewu/bugs.models/tree/master/nplcm);

- For Mac OS X system, one can install the package and study the components of
  each function. If you find package `ks` cannot be loaded due to failure of 
  loading package `rgl`, follow the following steps:
  
    1. install X11 by going [here](http://xquartz.macosforge.org/trac/wiki/X112.7.7)
    
    2. `install.packages("http://cran.r-project.org/src/contrib/rgl_0.95.1158.tar.gz",repo=NULL,type="source")`

