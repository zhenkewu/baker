#' baker: **B**ayesian **A**nalytic **K**it for 
#' **E**tiology **R**esearch 
#'
#' `baker` is designed for disease etiology studies from case-control data 
#' with multiple sources of measurements with potential errors. If you are 
#' interested in estimating the population etiology pie (a vector of fractions
#' that sum to one), and the probability of each cause for a particular 
#' individual case, try `baker`. 
#' 
#' `baker` implements hierarchical Bayesian models to infer disease etiology
#' for multivariate binary data. We created `baker` to catalyze effective 
#' communications between analysts and practicing clinicians that are vital to
#' the success of etiology studies. The `baker` package offers 
#' modules to  
#' \itemize{
#' \item Import and tidy the
#' PERCH data (the study that motivates the creation of this package), 
#' \item Transform, explore the data,
#' \item Specify, automatically generate the model files, and fit the models (npLCM),  
#' \item Store and visualize posterior summaries for communicating scientific 
#' findings, and
#' \item Check and compare the fitted models.
#' }
#' 
#' `baker` has implemented models for dependent 
#' measurements given disease status, regression analyses of etiology, 
#' multiple imperfect measurements, different priors for true positive rates 
#' among cases with differential measurement characteristics, and 
#' multiple-pathogen etiology. Scientists in Pneumonia Etiology Research for 
#' Child Health (PERCH) study usually refer to the etiology distribution 
#' as "population etiology pie" and "individual etiology pie" for their 
#' compositional nature, hence the name of the package (baking the pie).
#' 
#' 
#' @seealso 
#' \itemize{
#' \item <https://github.com/zhenkewu/baker> for the source code
#' and system/software requirements to use `baker` for your data.
#' }
#'
#' @import rjags
#' @return No returned value; documentation purpose only.
#' @section baker functions:
#' [nplcm()]
#'
#' @docType package
#' @name baker
NULL
#> NULL

