**baker**: Bayesian Analysis Kit for Etiology Research
------
> An R Package for Fitting Bayesian [Nested Partially Latent Class Models](http://biostats.bepress.com/jhubiostat/paper276/) 

[![Gitter](https://badges.gitter.im/Join Chat.svg)](https://gitter.im/zhenkewu/baker?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![Build Status](https://travis-ci.org/zhenkewu/baker.svg?branch=master)](https://travis-ci.org/zhenkewu/baker)

Installation
--------------
```r
install.packages("devtools",repos="https://cloud.r-project.org")
devtools::install_github("zhenkewu/baker")
```
Note: run `install.packages("pbkrtest")` for `R(>=3.2.3)` if this package is reported
as missing.

Vignettes
-------------
```r
devtools::install_github("zhenkewu/baker", build_vignettes=TRUE) # will take extra time to run a few examples.
browseVignettes("baker")
```

Graphical User Interface (GUI)
--------------------------------
```r
install.packages("devtools",repos="http://watson.nci.nih.gov/cran_mirror/")
devtools::install_github("zhenkewu/baker")
shiny::runGitHub("baker","zhenkewu",subdir="inst/shiny")
```

Analytic Goal
-------------------------------------

- To study disease etiology from case-control data from multiple sources that have measurement errors. If you are interested in estimating the population etiology pie (fraction), and the probability of each cause for individual case, try `baker`.


Comprison to Other Existing Solutions
------------------------------------------------
- Acknowledges various levels of measurement errors and combines multiple sources
of data for optimal disease diagnosis.

Main Functions
-----------------------------
- `nplcm()` that fits the model with or without covariates.


Details
-------------------------------------

1. Implements hierarchical Bayesian models to infer disease etiology for multivariate binary data. The package builds in functionalities for data cleaning, exploratory data analyses, model specification, model estimation, visualization and model diagnostics and comparisons, catalyzing vital effective communications between analysts and practicing clinicians. 
2. `baker` has implemented models for dependent measurements given disease status, regression analyses of etiology, multiple imperfect measurements, different priors for true positive rates among cases with differential measurement characteristics, and multiple-pathogen etiology.
3. Scientists in [Pneumonia Etiology Research for Child Health](http://www.jhsph.edu/research/centers-and-institutes/ivac/projects/perch/) (PERCH) study usually refer to the etiology distribution as "*population etiology pie*" and "*individual etiology pie*" for their compositional nature, hence the name of the package.
4. Reference publication can be found [here](http://onlinelibrary.wiley.com/doi/10.1111/rssc.12101/abstract) and [here](http://biostatistics.oxfordjournals.org/content/early/2016/08/19/biostatistics.kxw037.short?rss=1).

Platform
---------
- The `baker` package is compatible with OSX, Linux and Windows systems, each requiring a slightly different setup as described below. If you need to speed up the installation and analysis, please contact the 
maintainer or chat by clicking the `gitter` button at the top of this README file. 

Connect `R` to `JAGS` or `WinBUGS`
----------------------------------
#### Mac OSX 10.11 El Capitan

1. Use [Just Another Gibbs Sampler (JAGS)](http://mcmc-jags.sourceforge.net/)
2. Install JAGS 4.2.0; Download [here](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/)
3. Install `R`; Download from [here](https://cran.r-project.org/)
4. Fire up `R`, run `R` command `install.pacakges("rjags")`
5. Run `R` command `library(rjags)` in R console; If the installations are successfull, you'll see some notes like this:

    ```r
    >library(rjags)
    Loading required package: coda
    Linked to JAGS 4.x.0
    Loaded modules: basemod,bugs
    ```

- Run `R` command `library(baker)`. If the package `ks` cannot be loaded due to failure of loading package `rgl`, first install X11 by going [here](http://xquartz.macosforge.org/trac/wiki/X112.7.7), followed by
    ```r
    install.packages("http://download.r-forge.r-project.org/src/contrib/rgl_0.95.1504.tar.gz",repo=NULL,type="source")
    ```

#### Unix (Build from source without administrative privilege)

Here we use [JHPCE](https://jhpce.jhu.edu/) as an example. [The complete installation guide](https://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/) offers 
extra information. 

1. Download source code for [JAGS 4.2.0](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.2.0.tar.gz/download);

2. Suppose you've downloaded it in `~/local/jags/4.2.0`. Follow the bash commands below:

    ``` bash
    # decompress files:
    tar zxvf JAGS-4.2.0.tar.gz
    
    # change to the directory with newly decompressed files:
    cd ~/local/jags/4.2.0/JAGS-4.2.0
    
    # specify new JAGS home:
    export JAGS_HOME=$HOME/local/jags/4.2.0/usr
    export PATH=$JAGS_HOME/bin:$PATH
    
    # link to BLAS and LAPACK:
    # Here I have used "/usr/lib64/atlas/" and "/usr/lib64/" on JHPCE that give me
    # access to libblas.so.3 and liblapack.so.3. Please modify to paths on your system.
    LDFLAGS="-L/usr/lib64/atlas/ -L/usr/lib64/" ./configure --prefix=$JAGS_HOME --libdir=$JAGS_HOME/lib64 
    
    # if you have 8 cores:
    make -j8
    make install
    
    # prepare to install R package, rjags:
    export PKG_CONFIG_PATH=$HOME/local/jags/4.2.0/usr/lib64/pkgconfig 
    
    module load R
    R> install.packages("rjags")
    ```
    
3. Also check out the [INSTALLATION](https://cran.r-project.org/web/packages/rjags/INSTALL) file for `rjags` package.


#### Windows 

- JAGS 4.2.0
    1. Install `R`; Download from [here](https://cran.r-project.org/)
    2. Install [JAGS 4.2.0](https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/); Add the path to JAGS 4.2.0 into the environmental variable (essential for R to find the jags program). See [this](http://superuser.com/questions/949560/how-do-i-set-system-environment-variables-in-windows-10) for setting environmental variables;
    3. Fire up `R`, run `R` command `install.pacakges("rjags")`
    4. Install [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/) (for building and installing R pacakges from source); Add the path to `Rtools` (e.g. `C:\Rtools\`) into your environmental variables so that R knows where to find it. 

- [WinBUGS 1.4.3](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/)
    1. Install the [patch](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/the-bugs-project-winbugs-patches/)
    2. Install the WinBUGS 1.4.x [immortality key](http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/)

Maintainer:
--------------------------

Zhenke Wu (zhwu@jhu.edu)

Department of Biostatistics

Johns Hopkins Bloomberg School of Public Health
