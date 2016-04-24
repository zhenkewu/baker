**baker**: Bayesian Analysis Kit for Etiology Research
------
> An R Package for Fitting Bayesian [Nested Partially Latent Class Models](http://biostats.bepress.com/jhubiostat/paper276/) 

[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/zhenkewu/baker?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/zhenkewu/baker.svg?branch=master)](https://travis-ci.org/zhenkewu/baker)

How to install?
--------------
```r
install.packages("devtools",repos="http://watson.nci.nih.gov/cran_mirror/")
devtools::install_github("zhenkewu/baker")
```

How to run baker user interface?
--------------------------------
```r
install.packages("devtools",repos="http://watson.nci.nih.gov/cran_mirror/")
devtools::install_github("zhenkewu/baker")
shiny::runGitHub("baker","zhenkewu",subdir="inst/shiny")
```

Why should someone use `baker`?
-------------------------------------

- To study disease etiology from case-control data from multiple sources that have measurement errors. If you are interested in estimating the population etiology pie (fraction), and the probability of each cause for individual case, try `baker`.

Details
-------------------------------------

> * Implements hierarchical Bayesian models to infer disease etiology for multivariate binary data. The package builds in functionalities for data cleaning, exploratory data analyses, model specification, model estimation, visualization and model diagnostics and comparisons, catalyzing vital effective communications between analysts and practicing clinicians. 
  * `baker` has implemented models for dependent measurements given disease status, regression analyses of etiology, multiple imperfect measurements, different priors for true positive rates among cases with differential measurement characteristics, and multiple-pathogen etiology.
  * Scientists in [Pneumonia Etiology Research for Child Health](http://www.jhsph.edu/research/centers-and-institutes/ivac/projects/perch/) (PERCH) study usually refer to the etiology distribution as "*population etiology pie*" and "*individual etiology pie*" for their compositional nature, hence the name of the package.
    
- Reference publication can be found [here](http://onlinelibrary.wiley.com/doi/10.1111/rssc.12101/abstract) and [here](http://biostats.bepress.com/jhubiostat/paper276/).

How does it compare to other existing solutions?
------------------------------------------------
- Acknowledges various levels of measurement errors and combines multiple sources
of data.

What are the main functions?
-----------------------------
- `nplcm()` that fits the model with or without covariates.

Platform
---------
The `baker` package is compatible with OSX, Linux and Windows systems, each requiring a slightly different setup as described below. 

If you need to speed up the installation and analysis, please contact the 
maintainer or chat by clicking the `gitter` button at the top of this README file. 

Connect `R` to `JAGS` or `WinBUGS`
----------------------------------
- Mac OSX 10.11 El Capitan, or Linux on High Performance Computing facilities 
    - Use [Just Another Gibbs Sampler (JAGS)](http://mcmc-jags.sourceforge.net/)
    - Install JAGS 4.2.0; Download [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/)
    - Install `R`; Download from [here](https://cran.r-project.org/)
    - Fire up `R`, run `R` command `install.pacakges("rjags")`
    - Run `R` command `library(rjags)` in R console; If the installations are successfull, you'll see some notes like this:
    
    ~~~r
    >library(rjags)
    Loading required package: coda
    Linked to JAGS 4.x.0
    Loaded modules: basemod,bugs
    ~~~
    
    - Run `R` command `library(baker)`. If the package `ks` cannot be loaded due to failure of loading package `rgl`, follow the following steps:
          + install X11 by going [here](http://xquartz.macosforge.org/trac/wiki/X112.7.7);
          + `install.packages("http://download.r-forge.r-project.org/src/contrib/rgl_0.95.1200.tar.gz",repo=NULL,type="source")`
		  
- Windows (if using JAGS 4.2.0)
    + Install `R`; Download from [here](https://cran.r-project.org/)
    + Install [JAGS 4.2.0](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/); Add the path to JAGS 4.2.0 into the environmental variable (essential for R to find the jags program). See [this](http://superuser.com/questions/949560/how-do-i-set-system-environment-variables-in-windows-10) for setting environmental variables;
    + Fire up `R`, run `R` command `install.pacakges("rjags")`
	  + Install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/) (for building and installing R pacakges from source); Add the path to `Rtools` (e.g. `C:\Rtools\`) into your environmental variables so that R knows where to find it. 

- Windows (if using [WinBUGS 1.4.3](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/))
    + Install the [patch](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/the-bugs-project-winbugs-patches/)
    + Install the WinBUGS 1.4.x [immortality key](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/)

Maintainer:
--------------------------

Zhenke Wu (zhwu@jhu.edu)

Department of Biostatistics

Johns Hopkins Bloomberg School of Public Health
