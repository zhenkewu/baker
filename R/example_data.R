#' Hypothetical pathogens and their categories (virus or bacteria)
#'
#' This is used in simulations where the pathogen names are from the alphabet,
#' and we hope to plot etiologies grouped by virus or bacteria
#'
#' @format A matrix of two columns
#' \describe{
#'   \item{pathogen}{names of the hypothetical pathogens, A-Z}
#'   \item{pathogen_type}{category of the hypothetical pathogens, `B` for
#'   bacterium, `V` for virus, which are randomly assigned.}
#' }
#' @usage data("pathogen_category_simulation")
"pathogen_category_simulation"

#' Simulated dataset that is structured in the format necessary for an [nplcm()] with regression 
#'
#' Data set for illustrating regression functionalities
#' 
#' @format A list containing three items 
#' \describe{
#'   \item{Mobs}{BrS level measurements: N = 1,200 (half cases and half controls);
#'   one slice of BrS measurements (6 dimensional, A-F); one slice of SS measurements (2 dimensional,
#'   A and B)}
#'   \item{Y}{case-control status}
#'   \item{X}{matrix of covariates (N by 4); columns: SITE (1 and 2, each with 600 subjects), 
#'   DATE (index from 1:300), std_date (standardized DATE), ENRLDATE (actual dates)}
#' }
#' 
#' @usage data("data_nplcm_reg")
"data_nplcm_reg"


# NB: a few more data sets: data_nplcm_no_reg (no regression setting); this can be used in a lot of functions
# for visualization. Perhaps also need to add posterior samples as data set just as a way to 
# demonstrate the plotting functionality and unit test.